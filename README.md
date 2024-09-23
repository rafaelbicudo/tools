# tools
Collection of scripts for preparing, running and analyzing quantum chemistry simulations. Check the [tutorials](https://github.com/rafaelbicudo/tools/tree/main/tutorials) file for more information.

### fit_dihedral.py

Performs a linear regression to determine Ryckaert-Bellemans torsional coefficients used in the parameterization of classical dihedral energies. The optimization uses quantum mechanical calculations of a rigid scan via Gaussian09/16 code, performs a cubic interpolation in the quantum data and fit the coefficients to reproduce the polynomial. 

**_PS_**: As the [DICEtools](https://github.com/hmcezar/dicetools) package can only be installed via GitHub, the package path needs to be in the `$PATH` variable or the `fit_dihedral.py` script should be in the same directory of DICEtools scripts.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >= 3.10
* [numpy](https://numpy.org)
* [pandas](https://pandas.pydata.org)
* [sklearn](https://scikit-learn.org/stable/index.html)
* [scipy](https://scipy.org)
* [matplotlib](https://matplotlib.org)
* [dicetools](https://github.com/hmcezar/dicetools)

### Usage
```
$ python fit_dihedral.py -h
usage: fit_dihedral.py [-h] [--method METHOD] [--alpha ALPHA] [--weight WEIGHT] [--cutoff CUTOFF] [--remove-overlap] [--max-barrier MAX_BARRIER] [--fit-from-total]
                       [--angles ANGLES [ANGLES ...]]
                       gaussianlogfile xyzrotationsfile topfile txtfile dfrfile a1 a2 a3 a4 npoints

Performs a linear regression to determine torsional dihedral coefficients.

positional arguments:
  gaussianlogfile       the gaussian log file.
  xyzrotationsfile      file with all torsional scan configurations.
  topfile               the topology file.
  txtfile               DICE .txt file
  dfrfile               DICE .dfr file
  a1                    first atom defining the reference dihedral.
  a2                    second atom defining the reference dihedral.
  a3                    third atom defining the reference dihedral.
  a4                    fourth atom defining the reference dihedral.
  npoints               number of configurations during the scan.

options:
  -h, --help            show this help message and exit
  --method METHOD, -m METHOD
                        the method employed in the linear regression (least-square, ridge, lasso, lassocv).
  --alpha ALPHA         the coefficient multiplying L1 penalty in Lasso linear regression. Default is 0.01.
  --weight WEIGHT, -w WEIGHT
                        the weight given to total energy minima points. Default is 1.
  --cutoff CUTOFF, -c CUTOFF
                        minimum atomic distance tolerated. Default is 0.5.
  --remove-overlap, -r  remove the overlapping configurations.
  --max-barrier MAX_BARRIER, -b MAX_BARRIER
                        limit the torsional barriers to the provided value (in kcal/mol).
  --fit-from-total, -t  fit the torsional angle using the total energy.
  --angles ANGLES [ANGLES ...]
                        all dihedrals (in degrees) to be used in the fit (include the max and min angles).
```

## get_excited-states.py
Extract excited state properties from Gaussian09/16 output file and export them as a `.csv` file. 

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >=3.9

### Usage
```
$ python get_excited-states.py -h
usage: get_excited-states.py [-h] [-dir DIR] [-log LOGFILE [LOGFILE ...]] [-csv CSVFILE]

Extract configurations from Gromos87 (.gro) files.

options:
  -h, --help            show this help message and exit
  -dir DIR              Directory with Gaussian09/16 output files.
  -log LOGFILE [LOGFILE ...], --logfile LOGFILE [LOGFILE ...]
                        Gaussian09/16 output file.
  -csv CSVFILE, --csvfile CSVFILE
                        .csv data file name. Default is "exc_data.csv".
```

## prep_asec.py
Creates an input file from `.gro` configurations to perform [ASEC](https://doi.org/10.1016/j.cplett.2007.02.012) calculations using Gaussian09/16.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >=3.9

### Usage
```
$ python prep_asec.py -h
usage: prep_asec.py [-h] [--grofile GROFILE] [--configs_dir CONFIGS_DIR] [--atomnums ATOMNUMS [ATOMNUMS ...]] [--residues RESIDUES [RESIDUES ...]] [--resnums RESNUMS [RESNUMS ...]]
                    [--keywords KEYWORDS [KEYWORDS ...]] [--charge CHARGE] [--spin_multiplicity SPIN_MULTIPLICITY] [--output OUTPUT]
                    itpfile [itpfile ...]

Extract configurations from Gromos87 (.gro) files.

positional arguments:
  itpfile               topology file(s) (.itp) from GROMACS.

options:
  -h, --help            show this help message and exit
  --grofile GROFILE, -gro GROFILE
                        reference .gro configuration file.
  --configs_dir CONFIGS_DIR, -dir CONFIGS_DIR
                        path to the directory with .gro configurations.
  --atomnums ATOMNUMS [ATOMNUMS ...], -an ATOMNUMS [ATOMNUMS ...]
                        list of atom numbers treated with QM (e.g., 1-10 22 82-93).
  --residues RESIDUES [RESIDUES ...], -res RESIDUES [RESIDUES ...]
                        list of residues treated with QM.
  --resnums RESNUMS [RESNUMS ...], -rn RESNUMS [RESNUMS ...]
                        number of the residue(s) to be treated with QM. Default is 1.
  --keywords KEYWORDS [KEYWORDS ...], -k KEYWORDS [KEYWORDS ...]
                        Gaussian keywords for the calculation. Default is "B3LYP/6-31G(d,p) Charge".
  --charge CHARGE, -c CHARGE
                        total charge of the system. Default is 0.
  --spin_multiplicity SPIN_MULTIPLICITY, -ms SPIN_MULTIPLICITY
                        spin multiplicity of the system. Default is 1.
  --output OUTPUT, -o OUTPUT
                        name of the output file. Default is "asec.com".
```

## prep_qmmm.py
Create Gaussian09/16 input files to run calculations according to the s-QM/MM method. Molecules from the quantum region are saturated according to the [H link-atom approach](https://doi.org/10.1002/anie.200802019).

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >=3.10

### Usage
```
$ prep_qmmm.py [-h] [--grofile GROFILE] [--configs_dir CONFIGS_DIR] [--atomnums ATOMNUMS [ATOMNUMS ...]] [--residues RESIDUES [RESIDUES ...]] [--resnums RESNUMS [RESNUMS ...]]
                    [--link_scale_factor LINK_SCALE_FACTOR] [--n_neighbors N_NEIGHBORS] [--cutoff CUTOFF] [--keywords KEYWORDS [KEYWORDS ...]] [--charge CHARGE]
                    [--spin_multiplicity SPIN_MULTIPLICITY] [--output OUTPUT] [--test]
                    itpfile [itpfile ...]

Extract configurations from Gromos87 (.gro) files.

positional arguments:
  itpfile               topology file(s) (.itp) from GROMACS.

options:
  -h, --help            show this help message and exit
  --grofile GROFILE, -gro GROFILE
                        reference .gro configuration file.
  --configs_dir CONFIGS_DIR, -dir CONFIGS_DIR
                        path to the directory with .gro configurations.
  --atomnums ATOMNUMS [ATOMNUMS ...], -an ATOMNUMS [ATOMNUMS ...]
                        list of atom numbers treated with QM (e.g. 1-10 22 82-93).
  --residues RESIDUES [RESIDUES ...], -res RESIDUES [RESIDUES ...]
                        list of residues treated with QM.
  --resnums RESNUMS [RESNUMS ...], -rn RESNUMS [RESNUMS ...]
                        number of the residue(s) to be treated with QM. Default is 1.
  --link_scale_factor LINK_SCALE_FACTOR, -sf LINK_SCALE_FACTOR
                        link atom distance scale factor. Default is 0.71.
  --n_neighbors N_NEIGHBORS, -nn N_NEIGHBORS
                        number of closest neighbors to redistribute the charge. Default is 3.
  --cutoff CUTOFF, -cut CUTOFF
                        cutoff radius (in AA) to redistribute the charge.
  --keywords KEYWORDS [KEYWORDS ...], -k KEYWORDS [KEYWORDS ...]
                        Gaussian keywords for the calculation. Default is "B3LYP/6-31G(d,p) Charge".
  --charge CHARGE, -c CHARGE
                        total charge of the system. Default is 0.
  --spin_multiplicity SPIN_MULTIPLICITY, -ms SPIN_MULTIPLICITY
                        spin multiplicity of the system. Default is 1.
  --output OUTPUT, -o OUTPUT
                        name of the output file. Default is "calc.com".
  --test                If True, generates a .xyz file with partial charges as bismuth atoms for visualization.
```
