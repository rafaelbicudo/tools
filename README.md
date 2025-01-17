# tools
Collection of scripts for preparing, running and analyzing quantum chemistry simulations. Check [tutorials](https://github.com/rafaelbicudo/tools/blob/main/tutorials/) for more information:

Tutorial 1: [Fitting proper torsional constants for OPLS-AA force field](https://github.com/rafaelbicudo/tools/blob/main/tutorials/fit_dihedral.md)  
Tutorial 2: [Generating Gaussian input for QM/MM calculations](https://github.com/rafaelbicudo/tools/blob/main/tutorials/prep_qmmm.md)  

## fit_dihedral.py

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
usage: fit_dihedral2.py [-h] [--method METHOD] [--alpha ALPHA] [--weight WEIGHT] [--cutoff CUTOFF] [--remove-overlap]
                        [--max-barrier MAX_BARRIER] [--fit-from-total] [--angles ANGLES [ANGLES ...]]
                        [--etot_limits ETOT_LIMITS ETOT_LIMITS] [--edih_limits EDIH_LIMITS EDIH_LIMITS]
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
  --etot_limits ETOT_LIMITS ETOT_LIMITS
                        lower and upper limits of the total energy plot.
  --edih_limits EDIH_LIMITS EDIH_LIMITS
                        lower and upper limits of the dihedral energy plot.
```

## get_excited-states.py
Extract excited state properties from Gaussian09/16 output file and export them as a `.csv` file. 

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >=3.9

### Usage
```
$ python get_excited-states.py -h
usage: get_excited-states.py [-h] [-dir DIR] [-log LOGFILE [LOGFILE ...]] [-csv CSVFILE]

Extract excited state data from Gaussian output files.

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
usage: prep_asec.py [-h] [--grofile GROFILE] [--configs_dir CONFIGS_DIR] [--atomnums ATOMNUMS [ATOMNUMS ...]]
                    [--residues RESIDUES [RESIDUES ...]] [--resnums RESNUMS [RESNUMS ...]]
                    [--keywords KEYWORDS [KEYWORDS ...]] [--charge CHARGE] [--spin_multiplicity SPIN_MULTIPLICITY]
                    [--output OUTPUT]
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
$ python prep_qmmm.py -h
usage: prep_qmmm.py [-h] [--grofile GROFILE] [--configs_dir CONFIGS_DIR] [--atomnums ATOMNUMS [ATOMNUMS ...]]
                    [--residues RESIDUES [RESIDUES ...]] [--resnums RESNUMS [RESNUMS ...]]
                    [--link_scale_factor LINK_SCALE_FACTOR] [--n_neighbors N_NEIGHBORS] [--cutoff CUTOFF]
                    [--embedding_cutoff EMBEDDING_CUTOFF] [--percentage_redist PERCENTAGE_REDIST]
                    [--keywords KEYWORDS [KEYWORDS ...]] [--charge CHARGE] [--spin_multiplicity SPIN_MULTIPLICITY]
                    [--output OUTPUT] [--checkpoint] [--test]
                    itpfile [itpfile ...]

Generates input file for QM/MM calculations using Gaussian.

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
                        cutoff radius (in AA) to redistribute the charge from the removed atom.
  --embedding_cutoff EMBEDDING_CUTOFF, -ec EMBEDDING_CUTOFF
                        cutoff distance (in AA) to include point charges. Default is to not use it.
  --percentage_redist PERCENTAGE_REDIST, -pr PERCENTAGE_REDIST
                        percentage of distant point charges to redistribute the net charge from embedding cut. Default is 0.1
                        (i.e., 10%).
  --keywords KEYWORDS [KEYWORDS ...], -k KEYWORDS [KEYWORDS ...]
                        Gaussian keywords for the calculation. Default is "B3LYP/6-31G(d,p) Charge NoSymm".
  --charge CHARGE, -c CHARGE
                        total charge of the system. Default is 0.
  --spin_multiplicity SPIN_MULTIPLICITY, -ms SPIN_MULTIPLICITY
                        spin multiplicity of the system. Default is 1.
  --output OUTPUT, -o OUTPUT
                        name of the output file. Default is "calc.com".
  --checkpoint, -chk    If True, add a line to save the checkpoint file during calculations.
  --test                If True, generates a .xyz file for visualization, with partial charges as bismuth atoms.
```

## get_g16_opt_configs.py
Extract the transient configurations that were considered during a geometry optimization using Gaussian09/16.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >= 3.9

### Usage
```
$ python get_g16_opt_configs.py -h
usage: get_g16_opt_configs.py [-h] [-o OUTPUT] file

Get transient configurations from geometry optimization using Gaussian.

positional arguments:
  file                  Gaussian optimization output file.

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Name of the file with transient configurations. Default is opt_configs.xyz.
```

## make_interface.py
Duplicate a given molecule to create an interface. It was developed to create idealized aggregation types, _e.g._, H- and J-type.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >= 3.9
* [numpy](https://numpy.org)
* [MDAnalysis](https://www.mdanalysis.org)
* [sklearn](https://scikit-learn.org/stable/index.html)
* [scipy](https://scipy.org)

### Usage
```
$ python make_interface.py -h
usage: make_interface.py [-h] [-t TYPE] [-d DISTANCE] [-o OUTPUT] [-s SHEAR] [-rot ROT_ANGLE] file

Create interface between molecules.

positional arguments:
  file                  File with configuration parsed by MDAnalysis.

options:
  -h, --help            show this help message and exit
  -t TYPE, --type TYPE  Type of aggregation, i.e., H-, X-, T- or J-type. Default is H.
  -d DISTANCE, --distance DISTANCE
                        Distance between molecules. Default is 3 Angstrom.
  -o OUTPUT, --output OUTPUT
                        Name and format of the output file. Default is stacked_file.pdb.
  -s SHEAR, --shear SHEAR
                        Shear distance to shift the molecule for J-type aggregation. Default is 1.5 Angstrom.
  -rot ROT_ANGLE, --rot_angle ROT_ANGLE
                        Rotation angle for the X- or T-type aggregation. Default is 90 degrees.
```

## prep_glass_trans.py
Create input files for running GROMACS simulations at different temperatures to determine the glass transition temperature.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >= 3.9
* [numpy](https://numpy.org)

### Usage
```
$ python prep_glass_trans.py -h
usage: prep_glass_trans.py [-h] [-T0 INITIAL_TEMPERATURE] [-Tf FINAL_TEMPERATURE] [-dT TEMPERATURE_INTERVAL]
                           [-sub SUBMISSION_FILE]
                           gro_file top_file

Prepare GROMACS input files for the simulations required to determine the glass transition temperature.

positional arguments:
  gro_file              the .gro input file.
  top_file              the .top topology file.

options:
  -h, --help            show this help message and exit
  -T0 INITIAL_TEMPERATURE, --initial_temperature INITIAL_TEMPERATURE
                        the initial temperature (in Celsius). Default is 20.
  -Tf FINAL_TEMPERATURE, --final_temperature FINAL_TEMPERATURE
                        the final temperature (in Celsius). Default is 280.
  -dT TEMPERATURE_INTERVAL, --temperature_interval TEMPERATURE_INTERVAL
                        the temperature interval between simulations (in Celsius). Default is 5.
  -sub SUBMISSION_FILE, --submission_file SUBMISSION_FILE
                        file with the header of the submission file. Default is sub.txt.
```

## get_glass_trans.py
Determine the glass transition temperature from the output of the simulations prepared with `prep_glass_trans.py`.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >= 3.9
* [numpy](https://numpy.org)
* [pandas](https://pandas.pydata.org)
* [scipy](https://scipy.org)
* [matplotlib](https://matplotlib.org)

### Usage
```
$ python get_glass_trans.py -h
usage: get_glass_trans.py [-h] [-dirs DIRECTORIES [DIRECTORIES ...]] [-i VOLUME_INDEX] [-s] [-vd VOLUME_DATA]
                          [-pref DIRECTORY_PREFFIX] [-lT LOW_TEMPERATURES] [-hT HIGH_TEMPERATURES] [-dT TEMPERATURE_INTERVAL]
                          [-t_min MINIMUM_TEMPERATURE] [-t_max MAXIMUM_TEMPERATURE]

Determine the glass temperature transition from GROMACS output files.

options:
  -h, --help            show this help message and exit
  -dirs DIRECTORIES [DIRECTORIES ...], --directories DIRECTORIES [DIRECTORIES ...]
                        the directories for the independent sets of simulation. Used for more than 1 set of simulations. Default
                        is [].
  -i VOLUME_INDEX, --volume_index VOLUME_INDEX
                        the index for volume when running 'gmx energy'. Default is 22.
  -s, --save_data       If True, saves the volume data in the .csv format. Default is True.
  -vd VOLUME_DATA, --volume_data VOLUME_DATA
                        the volume data .csv file name (e.g., vol_data.csv). Default is None.
  -pref DIRECTORY_PREFFIX, --directory_preffix DIRECTORY_PREFFIX
                        the preffix of the temperature directories. Default is 'temp'.
  -lT LOW_TEMPERATURES, --low_temperatures LOW_TEMPERATURES
                        list with low temperature values to be included in the fitting. Default is 20-70.
  -hT HIGH_TEMPERATURES, --high_temperatures HIGH_TEMPERATURES
                        interval of high temperature values to be included in the fitting. Default is 230-280.
  -dT TEMPERATURE_INTERVAL, --temperature_interval TEMPERATURE_INTERVAL
                        the temperature interval between simulations (in Celsius). Default is 5.
  -t_min MINIMUM_TEMPERATURE, --minimum_temperature MINIMUM_TEMPERATURE
                        the minimum temperature for plotting (in Celsius). Default is 0.
  -t_max MAXIMUM_TEMPERATURE, --maximum_temperature MAXIMUM_TEMPERATURE
                        the maximum temperature for plotting (in Celsius). Default is 300.
```

## get_closest_mols.py
Extract snapshots from GROMACS simulations with a given atom group centered and the closest molecules to it, and export it to a .xyz file.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >= 3.9
* [numpy](https://numpy.org)
* [MDAnalysis](https://www.mdanalysis.org)

### Usage
```
$ python get_closest_mols.py -h
usage: get_closest_mols.py [-h] [-s TOPOLOGY_FILE] [-f TRAJECTORY_FILE] [-ref REFERENCE]
                           [-cref CLOSEST_MOLS_REFERENCE] [-n N_CLOSEST] [-init INITIAL_FRAME]
                           [-end FINAL_FRAME] [-n_frames NUMBER_OF_FRAMES] [-cm_to_origin] [-o OUTPUT_FILE]

Get closest molecules to a given reference from GROMACS output.

options:
  -h, --help            show this help message and exit
  -s TOPOLOGY_FILE, --topology_file TOPOLOGY_FILE
                        GROMACS structure+mass(db) file (.tpr, .gro, .g96, .pdb, .brk or .ent). Default is 'topol.tpr'.
  -f TRAJECTORY_FILE, --trajectory_file TRAJECTORY_FILE
                        GROMACS trajectory file (.xtc, .trr, .cpt, .gro, .g96, .pdb or .tng). Default is 'traj.xtc'.
  -ref REFERENCE, --reference REFERENCE
                        Atom selection language (same used in VMD) for the reference group. Default is 'resname UNK'.
  -cref CLOSEST_MOLS_REFERENCE, --closest_mols_reference CLOSEST_MOLS_REFERENCE
                        Atom selection language (same used in VMD) for the closest atoms. Default is 'not resname UNK'.
  -n N_CLOSEST, --n_closest N_CLOSEST
                        Number of closest molecules to be included. Default is 1.
  -init INITIAL_FRAME, --initial_frame INITIAL_FRAME
                        Initial frame. Default is 1.
  -end FINAL_FRAME, --final_frame FINAL_FRAME
                        Final frame. Default is -1 (last frame).
  -n_frames NUMBER_OF_FRAMES, --number_of_frames NUMBER_OF_FRAMES
                        Number of frames to be considered. Default is 100.
  -cm_to_origin, --center_of_mass_to_origin
                        If True, translates the center of mass of the reference molecule to the origin. Default is False.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file name. Default is 'output.xyz'.
```

### get_contacts.py

Compute the number of contacts between two AtomGroups from MDAnalysis.

### Dependencies
* [Python](https://scikit-learn.org/stable/index.html) >= 3.9
* [numpy](https://numpy.org)
* [pandas](https://pandas.pydata.org)
* [MDAnalysis](https://www.mdanalysis.org)

### Usage
```
$ python get_contacts.py -h
usage: get_contacts.py [-h] -g1 GROUP1 -g2 GROUP2 [-r RADIUS] [-o OUTPUT] topfile trajfile

Compute the number of contacts for molecules with backbone and side chains structures.

positional arguments:
  topfile               GROMACS topology/structure file path (e.g., .tpr or .gro).
  trajfile              GROMACS trajectory file path (e.g. .xtc or .trr).

options:
  -h, --help            show this help message and exit
  -g1 GROUP1, --group1 GROUP1
                        Atom selection language of Group 1.
  -g2 GROUP2, --group2 GROUP2
                        Atom selection language of Group 2.
  -r RADIUS, --radius RADIUS
                        Cutoff radius (in Angstrom) for contacts. Default is 4.5.
  -o OUTPUT, --output OUTPUT
                        Output file name. Default is 'contacts.csv'.
```