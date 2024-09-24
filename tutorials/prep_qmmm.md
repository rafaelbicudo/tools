# Generating Gaussian input for QM/MM calculations

This tutorial is dedicated to create Gaussian input files for QM/MM calculations with `prep_qmmm.py` script. The document is separated in two parts, for single and multiple files, respectively. 

## Single file

For this part we have 

In this example, we want the treat the central monomer as the quantum region and the rest as point charges.

## Multiple files

For this part we have 5 snapshots extracted from a Molecular Dynamics (MD) simulation of 5 oligomers in a solution of ethanol and water. The snapshots are located in the `configs` directory and the topology files are `ethanol.itp`, `tip3p.itp` and `ff_olig.itp`, such that:

```
$ ls
configs     ethanol.itp ff_olig.itp tip3p.itp
$ ls configs/
conf1_box1_000.gro conf1_box1_001.gro conf1_box1_002.gro conf1_box1_003.gro conf1_box1_004.gro
```

For example, the `conf1_box1_000.gro` is the following:

![5olig.tga](5olig.tga)

