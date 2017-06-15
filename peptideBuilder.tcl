#
# Peptide builder
#
# peptideBuilder.tcl,v 1.0 2017/04/25
#
#
#
# Distributed under GNU license
# Author:
#		Thomas Lemmin


package provide peptideBuilder 1.0


package require autopsf
package require utilities


namespace eval ::peptideBuilder:: {
	variable PBPATH [file dirname [file normalize [info script]]]
	variable aaa2a
	array set aaa2a {ALA A ARG R ASN N ASP D ASX B CYS C GLN Q GLU E
		GLX Z GLY G HIS H ILE I LEU L LYS K MET M PHE F PRO P SER S
		THR T TRP W TYR Y VAL V}
	variable a2aaa
	foreach rn [array names aaa2a] {
		set a2aaa($aaa2a($rn)) $rn
	}

	variable ideal_ss
	array set ideal_ss  {
			E {-120 113}
			G {-74 -4} 
			H {-57 -47}  
			T {-60 30}
			}

	proc ::peptideBuilder::build_peptide {fname file_type {name peptide}} {
		#This function builds coordinate for a sequence of amino acid (pdb and psf).
		#Parameters
		#	fname: path to file with sequence
		#	file_type:
		#		full
		#			format:
		#				resname phi psi
		#				Example:
		#				VAL -57 -47
		#				ALA -57 -47
		#		ideal
		#			format:
		#				ss,seq
		#				where:
		#			 		ss = {E, G, H, T}
		#					seq = sequence one letter code
		#				Example:
		#					H,ASDFWEFG
		#					E,SHTT
		#	name: name of pdb file
	
		variable PBPATH
		if {$file_type == "full"} {
			set atoms [build_full $fname]	
		} elseif {$file_type == "ideal"} {
			set atoms [build_ideal $fname]
		} else {
			puts "Unknown type: should be full or ideal"
			return -1
		}
		
		set numatoms [llength $atoms]
		set id [mol new atoms $numatoms]
		animate dup $id
		[atomselect $id all] set {name resname resid x y z} $atoms
		mol rename $id $name
		set topology [file join $PBPATH top_all36_prot.rtf]
		autopsf -mol $id -prefix $name -top $topology
		
		mol delete $id
	
	}
	
	proc ::peptideBuilder::build_full {fname} {
		set fp [open $fname r]
		set data [split [read $fp] "\n"]
		set atoms [list]
		set resid 1	
		set psi_previous 0
		foreach d $data {
			puts $d
			if {[llength $d] == 3} {
				lassign $d resname phi psi
				set atoms [build_backbone $atoms $phi $psi_previous $resname $resid]
				incr resid
				set psi_previous $psi
			} else {
				puts "Skipping line $d"
			}
		}
		return $atoms
	}
	
	proc ::peptideBuilder::build_ideal {fname} {
		variable a2aaa
		variable ideal_ss
	
		set fp [open $fname r]
		set data [split [read $fp] "\n"]
		set atoms [list]
		set resid 1	
		set ss_previous "L"
		foreach d $data {
			lassign [split $d ","] ss seq
			foreach s [split $seq ""] {
				set resname $a2aaa($s)
				if {[info exists ideal_ss($ss)]} {
					set phi [lindex $ideal_ss($ss) 0]
				} else {
					set phi 0
				}
			
				if {[info exists ideal_ss($ss_previous)]} {
					set psi [lindex $ideal_ss($ss_previous) 1]
				} else {
					set psi 0
				}
				puts "$ss $ss_previous $resname $phi $psi"
				set atoms [build_backbone $atoms $phi $psi $resname $resid]
				incr resid
				set ss_previous $ss
			}
			
		}

		return $atoms
	}
	proc ::peptideBuilder::build_backbone {atoms phi psi resname resid} {
		#atoms: name resname x y z
		#If atoms is empty, initialise 
		if {[llength $atoms] == 0} {
			set N_ {0 0 0}
			set CA_ {0 0 1}
			set C_ {0 1 1}
		} else {
			set N_ [lrange [lindex $atoms end-2] 3 end]
			set CA_ [lrange [lindex $atoms end-1] 3 end]
			set C_ [lrange [lindex $atoms end] 3 end]
		}
		foreach name {N CA C} {
			set atom [list]
			lappend atom $name
			lappend atom $resname
			lappend atom $resid
			if {$name == "N"} {
				set a1 $N_
				set a2 $CA_
				set a3 $C_
				set r  1.3558
				set theta 116.8400
				#psi
				set alpha $psi
			} elseif {$name == "CA"} {
				set a1 $CA_
				set a2 $C_
				set a3 [lrange [lindex $atoms end] 3 end]
				set r  1.4613
				set theta 126.7700
				set alpha	180	
			} elseif {$name == "C"} {
				set a1 $C_
				set a2 [lrange [lindex $atoms end-1] 3 end]
				set a3 [lrange [lindex $atoms end] 3 end]
				set r 1.5390
				set theta 114.4400
				#phi
				set alpha $phi
			} 
			set atom [concat $atom [build_cartesian $a1 $a2 $a3 $r $theta $alpha]]
			lappend atoms $atom
		}

		return $atoms
	}


	proc ::peptideBuilder::build_cartesian {a1 a2 a3 r theta phi} {
		#This function computes the cartesian coordinate of a missing atoms from provide internal coordinates
		#Parameters:
		#	a1: coordinates of first atom (list)
		#	a2: coordinates of second atom (list)	
		#	a3: coordinates of third atom (list)	
		#	r: coordinates of bond distance from a3 (Angstrom)	
		#	theta: angle between a2 a3 and missing atom (degrees)	
		#	phi: torsional angle between a1 a2 a3 and missing atom (degrees)	
		#Returns:
		#	newc: xyz coordinates of missing atom (list)

	  set theta [::util::deg2rad $theta] 
	  set phi   [::util::deg2rad $phi]
	  set cost  [expr cos($theta)]
	  set sint  [expr sin($theta)]
	  set cosp  [expr cos($phi)]
	  set sinp  [expr sin($phi)]
	  set rjk  [vecsub $a2  $a3]
	  set rjk [vecnorm $rjk]
	  set rij [vecsub $a1  $a2]

	  set cross [veccross $rij $rjk]
	  set cross [vecnorm $cross]
	  set cross2  [veccross $rjk $cross]
	  set cross2 [vecnorm $cross2]
	  set wt [list [expr $r * $cost] [expr $r * $sint * $cosp] [expr $r * $sint * $sinp]]
	  set newc [vecadd [vecscale $rjk [lindex $wt 0]] [vecscale $cross2 [lindex $wt 1]] [vecscale $cross [lindex $wt 2]]]
	  return [vecadd $a3 $newc]
	}
	
	proc ::peptideBuilder::extractPhiPsi {id {sel protein}} {
		set peptide [list]
		set a [atomselect $id "$sel and alpha"]
		set peptide [$a get {resname phi psi}]
		$a delete
		return $peptide
	}
}