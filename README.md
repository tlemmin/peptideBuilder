# peptideBuilder
This is a plug-in for VMD (requires autopsf). 
The code allows building a peptide ab initio just from a sequence and the corresponding secondary structure or backbone internal coordinates.

## Example
Input for an ideal secondary structure (e.g. helical peptide):
peptide.dat
```
h,ADFGEDT
```
Input for explicit torsion angles:
peptide_full.dat
```
ALA -57 -47
GLU -57 -47
PHE -57 -47
GLY -57 -47
```
Running it in VMD:
```
packae require peptideBuilder
::peptideBuilder::build_peptide peptide.dat
::peptideBuilder::build_full peptide_full.dat
```
