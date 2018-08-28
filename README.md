#X-raymol
---

## Description

This script written in Perl is a tool for identifying atom bonds, angles and dihedrals. 

## Usage

In order to execute this script:

```console
xraymol.pl input.xyz 
```

- input: a XYZ file the of the molecule.


## Typical output

For a 2-Me-pyrrole:

<pre>
------------------------------------
X-raymol
created by Evanildo Lacerda Jr.
------------------------------------

command -> $ xraymol [option] [file.xyz] 

[options]

none     generate .xray file
----     -------------------
-key     generate .key file [tinker]
-itp     generate .itp file [gromacs]
-h       help

# atoms = 13
# bonds = 13
# angles = 21
# torsions = 26
# bond criteria dij < 1.56

input: 2-Me-p_opt.xyz
output: 2-Me-p_opt.xray

</pre>



 
