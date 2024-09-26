#!/usr/bin/env python3

"""Creates an input file from .gro configurations to perform
   ASEC calculations using Gaussian09/16.
    
   Authors: Rafael Bicudo Ribeiro (@rafaelbicudo) and Leandro Rezende Franco
   DATE: 09/07/2024
"""

import argparse
import os
import sys


# Functions
def read_gro_file(grofile: str) -> dict:
    """Extract data from Gromos87 file.

    Adapted from "set_gaussian_qmmm.py"
    (https://github.com/rafaelbicudo/PhD).

    Args:
        grofile (str): .gro file name.

    Returns:
        data (dict): data from .gro file. 
    """
    
    data = {"resnum"   : [],
            "resname"  : [],
            "atomname" : [],
            "atomnum"  : [],
            "x_coord"  : [],
            "y_coord"  : [],
            "z_coord"  : []
            }

    # Read the grofile
    with open(grofile, "r") as f:
        lines = f.readlines()

        # Get the header
        data["header"] = lines[0]

        # Get the number of atoms
        data["n_atoms"] = int(lines[1].split()[0])
        
        # Get the box dimensions
        data["box_dim"] = (float(lines[-1].split()[0]), 
                           float(lines[-1].split()[0]), 
                           float(lines[-1].split()[0]))

        for line in lines[2:-1]:

            # Get the residue number
            data["resnum"].append(int(line[:5]))

            # Get the residue name
            data["resname"].append(line[5:10].split()[0])

            # Get the atom name
            data["atomname"].append(line[10:15].split()[0])

            # Get the atom number
            data["atomnum"].append(int(line[15:20]))

            # Get the atomic coordinates
            data["x_coord"].append(float(line[20:28]))
            data["y_coord"].append(float(line[28:36]))
            data["z_coord"].append(float(line[36:44]))

    return data


def read_itp_file(itpfile: str) -> dict:
    """Extract data from GROMACS topology (.itp) file.

    Adapted from "set_gaussian_qmmm.py"
    (https://github.com/rafaelbicudo/PhD).

    Args:
        itpfile (str): .itp file name.

    Returns:
        data (dict): data from .itp file.
    """

    def parse_file(f) -> None:
        """Reads .itp files and parse the data.

        Args:
            f (_io.TextIOWrapper): opened file to be read.
        """
        line = f.readline()

        # Get the atomic data
        while "[ atoms ]" not in line: 
            line = f.readline()

        while line:
            line = f.readline()
            if line.strip() == "" or line.startswith("[ "):
                break 
            elif line.startswith(";"):
                pass
            else:
                words = line.split()

                # Get the atom number
                data["atomnum"].append(words[0])

                # Get the atom type
                data["atomtype"].append(words[1])

                # Get the residue name
                data["resname"].append(words[3])

                # Get the atom name
                data["atomname"].append(words[4])

                # Get the atomic charge    
                data["atomcharge"].append(float(words[6]))

    data = {"atomnum"    : [],
            "atomtype"   : [],
            "resname"    : [],
            "atomname"   : [],
            "atomcharge" : []
            }

    # Read the topology file(s)
    if isinstance(itpfile, list):
        for file in itpfile:
            with open(file, "r") as f:
                parse_file(f)

    else:
        with open(itpfile, "r") as f:
            parse_file(f)
                
    return data


def get_QM_atoms(
    grofile: str, 
    atomnums: list,
    residues: list,
    resnums: int
) -> list:
    """Get a list with atoms to be treated with QM.

    Adapted from "get_configs_from_gro.py "
    (https://github.com/rafaelbicudo/PhD).

    Args:
        grofile (str): .gro file name.
        atomnums (list): QM atom numbers (indexes) from the .gro file.
        residues (list): QM residues from the .gro file.
        resnums (int): number of the residue to get the QM atoms from .gro file.
                       Only used when a list of residues is provided.

    Returns:
        qm_atoms (list): list with QM atoms. 
    """

    # Get the data dictionary
    gro_data = read_gro_file(grofile)

    # Checking for each combination of input to create a list with QM atoms
    if atomnums and not residues:
        qm_atoms = atomnums

        # Get the atom type and residue of each atom
        res = []
        for i in atomnums:
            for j in range(len(gro_data["resname"])):
                if gro_data["atomnum"][j] == i:
                    res.append(gro_data["resname"][j])

    elif not atomnums and residues:
        print(f"Selecting residue names {residues} with residue numbers {resnums}.\n")

        qm_atoms = []
        res = []
        for i in residues:
            for j in range(len(gro_data["resname"])):
                if gro_data["resname"][j] == i:
                    for k in resnums:
                        if k == gro_data["resnum"][j]:
                            qm_atoms.append(gro_data["atomnum"][j])
                            res.append(i)

    else:
        print('Please provide atom numbers or residues for the QM atoms.')
        sys.exit()

    return qm_atoms, res


def get_charge(
    itpfile: list,
    atom_name: str,
    res_name: str
) -> float:
    """Get the charge from the .itp file(s).

    Adapted from "get_configs_from_gro.py "
    (https://github.com/rafaelbicudo/PhD).

    Args:
        itpfile (str): list with .itp file names.
        atom_name (str): name of the atom.
        res_name (str): name of the residue.

    Returns:
        charge (float): charge of "atom_name".
    """

    # Loop over all .itp files
    for file in itpfile:

        # Get the data dictionary
        itp_data = read_itp_file(file)

        # Loop over all atom names
        for i in range(len(itp_data["atomname"])):
            
            # Check for matching atom names
            if itp_data["atomname"][i] == atom_name:

                # Check for matching residue names
                if itp_data["resname"][i] == res_name:

                    charge = itp_data["atomcharge"][i]
                    return charge
                
    print("Couldn't find the atom in the topology file(s).")
    sys.exit()


def write_asec_format(
    itpfile: list[str],
    qm_atoms: list,
    configs_dir: str,
    keywords: str,
    charge: int,
    spin_mult: int,
    output: str
) -> None:
    """Write the file for running ASEC on Gaussian09/16.

    Args:
        itpfile (list[str]): .itp file name.
        qm_atoms (list): list with QM atoms.
        configs_dir (str): path to the directory with .gro configurations.
        keywords (str): Gaussian keywords for the calculation.
        charge (int): total charge of the system.
        spin_mult (int): spin multiplicity of the system.
        output (str, optional): name of the output file.

    Returns:
        A Gaussian09/16 input file.
    """

    # Get the number of configurations
    n_configs = 0
    for file in os.listdir(configs_dir):
        if file.endswith(".gro"):
            n_configs += 1

    # Print the amount of configurations
    print(f"Considering {n_configs} configurations found in '{configs_dir}' directory.\n")

    # Write the coordinates for the Gaussian input file
    with open(f"{output}", "w") as fout:

        # Write the header
        fout.write(f"chk={output.split('.')[0]}.chk \n")
        fout.write(f"#p {' '.join(keywords)} \n\n")
        fout.write(f"ASEC calculation with {n_configs} configurations from '{configs_dir}' directory \n\n")
        fout.write(f"{charge} {spin_mult}\n")

        # Loop over all configuration files
        first = True
        for f in os.listdir(configs_dir):
            # Track the current configuration
            print(f"Working on configuration: {f}")

            # Get the data dictionaries
            gro_data = read_gro_file(os.path.join(configs_dir, f))

            # Write the QM coordinates
            if first:
                for j in range(len(gro_data["resname"])):
                    if gro_data["atomnum"][j] in qm_atoms:
                        
                        # Check if the second letter of atom name is lower case
                        if len(gro_data["atomname"][j])>1 and gro_data["atomname"][j][1].islower():
                            atom_name = gro_data["atomname"][j][:2]
                        else:
                            atom_name = gro_data["atomname"][j][:1]
                        
                        # Write the atom type and coordinates
                        fout.write("{}\t{:>.3f}\t{:>.3f}\t{:>.3f}\n".format(
                                    atom_name, gro_data["x_coord"][j]*10, 
                                    gro_data["y_coord"][j]*10, gro_data["z_coord"][j]*10))
                        
                fout.write("\n")

                # Write the other atoms as point charges
                for j in range(len(gro_data["resname"])):
                    if gro_data["atomnum"][j] not in qm_atoms:

                        # Get the charge
                        atom_charge = get_charge(itpfile, gro_data["atomname"][j], gro_data["resname"][j])

                        # Write the coordinates and partial charge
                        fout.write("{:>.3f}\t{:>.3f}\t{:>.3f}\t{:>.4f}\n".format(
                                    gro_data["x_coord"][j]*10, gro_data["y_coord"][j]*10,
                                    gro_data["z_coord"][j]*10, atom_charge/n_configs))

                first = False

            else:
                # Write the other atoms from other configurations as point charges
                for j in range(len(gro_data["resname"])):
                    if gro_data["atomnum"][j] not in qm_atoms:

                        # Get the charge
                        atom_charge = get_charge(itpfile, gro_data["atomname"][j], gro_data["resname"][j])

                        # Write the coordinates and partial charge
                        fout.write("{:>.3f}\t{:>.3f}\t{:>.3f}\t{:>.4f}\n".format(
                                    gro_data["x_coord"][j]*10, gro_data["y_coord"][j]*10,
                                    gro_data["z_coord"][j]*10, atom_charge/n_configs))
            
            fout.write("\n")

    # Print a final message
    print(f"\nGaussian input file \'{output}\' was successfully created.")

def parse_range(value: str):
    """Combine values from nested lists into a single list.

    Args:
        value (str): string with intervals

    Returns:
        (list): list with all values in a single list.
    """
    # Check if the input is in the form of a range (e.g., '1-100')
    if '-' in value:
        start, end = value.split('-')
        return list(range(int(start), int(end) + 1))
    else:
        # Otherwise, it might be a single integer
        return [int(value)]


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract configurations from Gromos87 (.gro) files.")

    parser.add_argument("itpfile", help="topology file(s) (.itp) from GROMACS.", nargs="+", type=str, default=[])
    parser.add_argument("--grofile", "-gro", help="reference .gro configuration file.", type=str)
    parser.add_argument("--configs_dir", "-dir", help="path to the directory with .gro configurations.", type=str)
    parser.add_argument("--atomnums", "-an", help="list of atom numbers treated with QM (e.g., 1-10 22 82-93).", nargs="+", type=parse_range, default=[])
    parser.add_argument("--residues", "-res", help="list of residues treated with QM.", nargs="+", type=str, default=[])
    parser.add_argument("--resnums", "-rn", help="number of the residue(s) to be treated with QM. Default is 1.", nargs="+", type=int, default=[1])
    parser.add_argument("--keywords", "-k", help="Gaussian keywords for the calculation. Default is \"B3LYP/6-31G(d,p) Charge\".", nargs="+", type=str, default=["B3LYP/6-31G(d,p) Charge"])
    parser.add_argument("--charge", "-c", help="total charge of the system. Default is 0.", type=int, default=0)
    parser.add_argument("--spin_multiplicity", "-ms", help="spin multiplicity of the system. Default is 1.", type=int, default=1)
    parser.add_argument("--output", "-o", help="name of the output file. Default is \"asec.com\".", type=str, default="asec.com")

    args = parser.parse_args()

    if args.grofile and args.configs_dir:
        print("Choose between a single '-gro' file or the directory 'configs_dir' with more than one file.")
        sys.exit()
    elif not args.grofile and not args.configs_dir:
        print("Provide a '-gro' file or a 'configs_dir' directory with files.")
        sys.exit()

    # Creates a single list of with atom numbers, to enable range formats, e.g. "1-10 25 33-38".
    allatomnums = []
    if args.atomnums:
        for sublist in args.atomnums:
            allatomnums.extend(sublist)

    if args.configs_dir:
        # Get the list of atoms treated with QM
        for root, _, files in os.walk(args.configs_dir):
            # Use the first file as reference file
            ref_file = os.path.join(root, files[0])

            qm_atoms, _ = get_QM_atoms(ref_file,
                                       allatomnums,
                                       args.residues,
                                       args.resnums)
    else:
        # Get the list of atoms treated with QM
        qm_atoms, _ = get_QM_atoms(args.grofile, allatomnums, args.residues, args.resnums)

    # Write the Gaussian input file
    write_asec_format(args.itpfile, qm_atoms, args.configs_dir, args.keywords, args.charge, args.spin_multiplicity, args.output)


if __name__ == '__main__':
    main()