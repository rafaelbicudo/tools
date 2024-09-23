#!/usr/bin/env python3

"""Create Gaussian09/16 input files to run calculations according to the s-QM/MM method.

    [X] Add an option to loop over a directory with several files.
    [ ] Add a flag to save the .chk file.
    [X] Fix the error when trying to run '-dir' in the cluster.
    [ ] Fix an issue of --test that creates a ".xyz" file if the a termination (e.g., .com) is not specified.

    AUTHOR: Rafael Bicudo Ribeiro (@rafaelbicudo)
    DATE: 08/2024
"""


import argparse
import sys
import os


# Functions
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


def read_gro_file(grofile: str) -> dict:
    """Extract data from Gromos87 file.

    Args:
        grofile (str): .gro file name.

    Returns:
        data (dict): data from .gro file. 
    """
    
    data = {"resnum"     : [],
            "resname"    : [],
            "atomname"   : [],
            "atomnum"    : [],
            "x_coord"    : [],
            "y_coord"    : [],
            "z_coord"    : [],
            "itpAtomNum" : []
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

    # Get the itp ordering by tracking the change in residue number
    for i in range(len(data["resnum"])):
        if (i == 0 or data["resnum"][i] != data["resnum"][i-1]):
            data["itpAtomNum"].append(1)
            j = 2
        else:
            data["itpAtomNum"].append(j)
            j += 1

    return data


def read_itp_file(itpfile: str) -> dict:
    """Extract data from GROMACS topology (.itp) file.

    Args:
        itpfile (str): .itp file name.

    Returns:
        data (dict): data from .itp file.
    """

    def parse_file(f) -> None:
        """Reads .itp files and parse the data.

        Args:
            f (_io.TextIOWrapper): _description_
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
                data["atomnum"].append(int(words[0]))

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

    Args:
        grofile (str): .gro file name.
        itpfile (list[str]): list with .itp file names.
        atomnums (list): QM atom numbers (indexes) from the .gro file.
        residues (list): QM residues from the .gro file.
        resnums (int): number of the residue to get the QM atoms from .gro file.
                       Only used when a list of residues is provided.

    Returns:
        qm_atoms (list): list with QM atoms with atom numbers from .gro file. 
    """

    # Get the data dictionary
    gro_data = read_gro_file(grofile)

    # Checking for each combination of input to create a list with QM atoms
    if atomnums and not residues:
        qm_atoms = atomnums

    elif not atomnums and residues:
        print(f"Selecting residue names {residues} with residue numbers {resnums}.\n")

        qm_atoms = []
        for i in residues:
            for j in range(len(gro_data["resname"])):
                if gro_data["resname"][j] == i:
                    for k in resnums:
                        if k == gro_data["resnum"][j]:
                            qm_atoms.append(gro_data["atomnum"][j])

    else:
        print('Please provide atom numbers (from .gro file) or residues of the QM atoms.')
        sys.exit()

    return qm_atoms


def get_charge(
    itpfile: str,
    atom_name: str,
    res_name: str,
    itpAtomNum: int
) -> float:
    """Get the charge from the .itp file(s).

    Args:
        itpfile (str): .itp file name.
        atom_name (str): name of the atom.
        res_name (str): name of the residue.
        itpAtomNum (int): corresponding atom number in the .itp file.
        
    Returns:
        charge (float): charge of "atom_name".
    """

    # Get the data dictionary
    itp_data = read_itp_file(itpfile)

    # Loop over all atom names
    for i in range(len(itp_data["atomname"])):
        # Check for matching itpAtomNumbers
        if itp_data["atomnum"][i] == itpAtomNum:

            # Check for matching atom names
            if itp_data["atomname"][i] == atom_name:

                # Check for matching residue names
                if itp_data["resname"][i] == res_name:

                    charge = itp_data["atomcharge"][i]
                    # print(f"Charge of {charge} found in the {itpfile} topology file.")

                    return charge
    
    # print(f"Couldn't find the atom in the {itpfile} topology file.")
    # sys.exit()


def get_connectivity(
    itpfile: list, 
    atomnum: int, 
    resname: str
) -> list:
    """Returns a list with bonded atoms to atomnum.

    Args:
        itpfile (list[str]): list with .itp file names.
        atomnum (int): atom number in the .itp file.
        resname (str): residue name of the corresponding atom.

    Returns:
        bonded_atoms (list): list with bonded atoms numbers (from the .itp file)
                             to the corresponding atom.
    """

    # Start the bonded atoms numbers list
    bonded_atoms = []

    # Open the .itp file(s)
    for file in itpfile:
        
        # Get the .itp data
        itp_data = read_itp_file(file)

        # Check if the residue is in the current .itp file
        res = False
        if resname in set(itp_data["resname"]):
            res = True
        
        # If not, save time by going to the next file
        if not res:
            continue

        with open(file, "r") as f:
            line = f.readline()

            # Read until the bonds block is found
            while "[ bonds ]" not in line:
                line = f.readline()

            # Reads until the file or the bonds block end
            while line:
                line = f.readline()
                if line.strip() == "" or line.startswith("[ "):
                    break 

                words = line.split()
                if str(f"{atomnum}") == words[0]: 
                    bonded_atoms.append(words[1])
                elif str(f"{atomnum}") == words[1]:
                    bonded_atoms.append(words[0])

        # Leaves the loop if the residue was found
        if res:
            break

    # Convert the bonded atom numbers into integers
    bonded_atoms = [int(x) for x in bonded_atoms]

    return file, bonded_atoms


def get_closest_atoms(
    grofile: str,
    target: int,
    cutoff: float,
    n_atoms: int
) -> list:
    """Returns a list of the n_atoms-th closest atoms.

    Args:
        grofile (str): .gro file with all atoms.
        target (int): atom number of the target atom to get.
        cutoff (float): cutoff distance (in AA).
        n_atoms (int): number of closest atoms to be returned. 

    Returns:
        closest_atoms (list): list with atom numbers of the closest atoms.
    """

    # Get the data dictionary
    gro_data = read_gro_file(grofile)

    # Define the reference coordinates
    x = gro_data["x_coord"][gro_data["atomnum"].index(target)]
    y = gro_data["y_coord"][gro_data["atomnum"].index(target)]
    z = gro_data["z_coord"][gro_data["atomnum"].index(target)]

    # Compute the distance with all atoms from .gro file
    distances = []
    for x_, y_, z_ in zip(gro_data["x_coord"], gro_data["y_coord"], gro_data["z_coord"]):
        # Compute the euclidean distance and append it to the list
        dist = ((x_-x)**2+(y_-y)**2+(z_-z)**2)**(1/2)
        distances.append(dist)

    # Create a list with atomic numbers
    atom_nums = gro_data["atomnum"].copy()

    # Sort the atomic numbers with increasing distances
    sorted_list = sorted(list(zip(distances, atom_nums)), key=lambda x: x[0])

    # Unzip the paired lists
    distances, closest_atoms = zip(*sorted_list)

    # Convert back to list
    distances = list(distances)
    closest_atoms = list(closest_atoms)

    if cutoff == 0.0:
        # Return the first "n_atoms" from this list
        return closest_atoms[1:n_atoms+1], distances[1:n_atoms+1]
    else:
        # Return the atoms closer than the "cutoff" distance
        for i in range(len(distances)):
            if distances[i] > cutoff/10:
                return closest_atoms[1:i], distances[1:i]


def h_link_saturation(
    grofile: str,
    itpfile: str,
    qm_atoms: list,
    dist_scale_factor: float
) -> None:
    """Saturates the QM region with hydrogens.

    Args:
        grofile (str): .gro file with all atoms.
        itpfile (str): .itp topology file(s).
        qm_atoms (list): list with atoms to be treated with QM.
        dist_scale_factor (float): scale factor to the link atom bond distance.

    Returns:
        to_remove_num (list): list with .gro atom numbers to be removed.
        h_links (list): list with coordinates of hydrogen link atoms.
    """

    # Get the data dictionaries
    gro_data = read_gro_file(grofile)

    # Initialize the list of atoms to be removed
    to_remove_num = []
    link_coords = []

    # Loop over QM atoms
    for i in qm_atoms:

        # Get the resname and the atom number (in the .itp file) of each QM atom
        itpAtomNum = gro_data["itpAtomNum"][gro_data["atomnum"].index(i)]
        resname = gro_data["resname"][gro_data["atomnum"].index(i)]
        resnum = gro_data["resnum"][gro_data["atomnum"].index(i)]

        # Get the bonded atoms (.itp file's atom number) to the QM atom
        _, bonded_atoms = get_connectivity(itpfile, itpAtomNum, resname)

        # Loop over the 20 closest atoms 
        closest_atoms, _ = get_closest_atoms(grofile, i, cutoff=0, n_atoms=20)

        for atom in closest_atoms:
            # Get the resname and the atom number (in the .itp file) of the neighbor
            _itpAtomNum = gro_data["itpAtomNum"][gro_data["atomnum"].index(atom)]
            _resnum = gro_data["resnum"][gro_data["atomnum"].index(atom)]

            # Check if the neighbor is bonded to the QM atom and 
            # belongs to different residues to find the atoms to be removed
            if (_itpAtomNum in bonded_atoms
                and _resnum == resnum
                and atom not in qm_atoms):
                # Add the atom to the removal list
                to_remove_num.append(atom)

                # Get the atomic coordinates of both QM and classical atoms
                x = gro_data["x_coord"][gro_data["atomnum"].index(i)]
                y = gro_data["y_coord"][gro_data["atomnum"].index(i)]
                z = gro_data["z_coord"][gro_data["atomnum"].index(i)]

                x_ = gro_data["x_coord"][gro_data["atomnum"].index(atom)]
                y_ = gro_data["y_coord"][gro_data["atomnum"].index(atom)]
                z_ = gro_data["z_coord"][gro_data["atomnum"].index(atom)]

                # Determine the link atom coordinates
                x_l = x + dist_scale_factor * (x_ - x)
                y_l = y + dist_scale_factor * (y_ - y)
                z_l = z + dist_scale_factor * (z_ - z)

                link_coords.append((x_l, y_l, z_l))

    return to_remove_num, link_coords


def get_charge_shifts(
    grofile: str,
    itpfile: str,
    qm_atoms: list,
    to_remove: list,
    n_neighbors: int,
    cutoff: float
) -> float | list:
    """Compute the charges to neutralize the system.

    Args:
        grofile (str): .gro file with all atoms.
        itpfile (str): .itp topology file(s).
        qm_atoms (list): list with atoms to be treated with QM.
        to_remove (list): list with atoms to be removed.
        n_neighbors (int): number of closest neighbors to redistribute the charge.
        cutoff (float): cutoff radius to redistribute the charge.
    Returns:
        qm_charge_shift (float): charge shift per atom from the QM atoms.
        neigh_sums (list): list with neighbor atoms.
        neigh_charge_shifts (list): list with charges to be added to the neighbor atoms.
    """

    # Initialize the variables
    qm_charge = 0.0
    neigh_nums = []
    neigh_charge_shifts = []

    # Get the data dictionary
    gro_data = read_gro_file(grofile)

    # Get the total charge of the QM atoms
    for j in qm_atoms:
        # Get the atom and residue name of each to be removed
        atomname = gro_data["atomname"][gro_data["atomnum"].index(j)]
        resname = gro_data["resname"][gro_data["atomnum"].index(j)]
        itpAtomNum = gro_data["itpAtomNum"][gro_data["atomnum"].index(j)]

        # Loop over the .itp files
        for file in itpfile:
            charge = get_charge(file, atomname, resname, itpAtomNum)
            if charge:
                qm_charge += charge 
                break

    # Loop over the atoms to be removed
    for i in to_remove:
        # Get the names and .itp number of each atom to be removed
        atomname = gro_data["atomname"][gro_data["atomnum"].index(i)]
        resname = gro_data["resname"][gro_data["atomnum"].index(i)]
        itpAtomNum = gro_data["itpAtomNum"][gro_data["atomnum"].index(i)]

        # Gets the charge of the atom
        for file in itpfile:
            # Search for the charge
            charge = get_charge(file, atomname, resname, itpAtomNum)

            # Breaks the loop if the charge is found
            if charge:
                break

        # Get all the atoms inside a sphere with "cutoff" radius
        closest_atoms, distances = get_closest_atoms(grofile, i, cutoff, n_neighbors)

        # Get the the "n_neighbors" closest atoms
        closest_atoms = closest_atoms[:n_neighbors]
        distances = distances[:n_neighbors]

        # Add the closest atoms to the list with neighbors
        neigh_nums.extend(closest_atoms)

        # Equally redistribute the charge over the "n_neighbors" atoms
        for i in closest_atoms:
            neigh_charge_shifts.append(charge/len(closest_atoms))

    # Combine charge shifts for repeated atoms
    summed_charges = {}
    for i, j in zip(neigh_nums, neigh_charge_shifts):
        if i in summed_charges:
            summed_charges[i] += j
        else:
            summed_charges[i] = j

    neigh_nums = list(summed_charges.keys())
    neigh_charge_shifts = list(summed_charges.values())

    # Compute the overall charge shift
    qm_charge_shift = qm_charge/(gro_data["n_atoms"]-len(to_remove)-len(qm_atoms))

    return qm_charge_shift, neigh_nums, neigh_charge_shifts


def write_gaussian_input(
    grofile: str,
    itpfile: str,
    qm_atoms: list,
    to_remove: list,
    link_coords: list,
    n_neighbors: int,
    cutoff: float,
    keywords: str,
    charge: int,
    spin_mult: int,
    output: str,
    test: bool,
) -> None:
    """Write the Gaussian09/16 input file.

    Args:
        grofile (str): .gro file with all atoms.
        itpfile (str): .itp topology file(s).
        qm_atoms (list): list with atoms to be treated with QM.
        to_remove (list): list with atoms to be removed.
        link_coords (list): list with link atoms coordinates.
        n_neighbors (int): number of closest neighbors to redistribute the charge.
        cutoff (float): cutoff radius to redistribute the charge.
        keywords (str): calculation keywords (e.g. HF/STO-3G Charge)
        charge (int): system's total charge.
        spin_mult (int): system's spin multiplicity.
        output (str): name of the output file.
        test (bool): write point charges as bismuth atoms for visualization.
    """

    # Get the data dictionaries
    gro_data = read_gro_file(grofile)

    # Get the charge shifts
    qm_shift, neigh_nums, neigh_shifts = get_charge_shifts(
        grofile, 
        itpfile, 
        qm_atoms, 
        to_remove, 
        n_neighbors, 
        cutoff
    )

    # Write the coordinates in the Gaussian input file
    with open(f"{output}", "w") as fout:
        if test:
            # Write the amount of atoms for the .xyz file
            fout.write(f"{len(gro_data['atomnum'])}\n\n")

        else:
            # Write the header
            fout.write(f"%chk={output.split('.')[0]}.chk \n")
            fout.write(f"#p {' '.join(keywords)} \n\n")
            fout.write(f"QM calculation with point charges \n\n")
            fout.write(f"{charge} {spin_mult}\n")

        # Write the QM atoms
        for j in range(len(gro_data["resname"])):
            # Write the QM atoms
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
        
        # Write the QM link atoms
        for j in link_coords:
            fout.write("H\t{:>.3f}\t{:>.3f}\t{:>.3f}\n".format(j[0]*10, j[1]*10, j[2]*10))

        # Write the blank line required by Gaussian
        if not test:
            fout.write("\n")
            
        # Write the remaining atoms as point charges
        for j in range(len(gro_data["resname"])):
            # Write the other atoms as point charges
            if (gro_data["atomnum"][j] not in qm_atoms
                and gro_data["atomnum"][j] not in to_remove):

                # Get the charge
                _charge = get_charge(itpfile, gro_data["atomname"][j], gro_data["resname"][j], gro_data["itpAtomNum"][j])

                # Redistribute the charges to neutralize the system
                if gro_data["atomnum"][j] in neigh_nums:
                    _charge += qm_shift + neigh_shifts[neigh_nums.index(gro_data["atomnum"][j])]
                else:
                    _charge += qm_shift

                # Boolean variable to check if partial charges are correctly placed
                if test:
                    # Write the partial charges as bismuth atoms
                    fout.write("Bi\t{:>.3f}\t{:>.3f}\t{:>3f}\n".format(
                                gro_data["x_coord"][j]*10, 
                                gro_data["y_coord"][j]*10,
                                gro_data["z_coord"][j]*10))

                else:
                    # Write the coordinates and partial charge
                    fout.write("{:>.3f}\t{:>.3f}\t{:>.3f}\t{:>}\n".format(
                                gro_data["x_coord"][j]*10, 
                                gro_data["y_coord"][j]*10,
                                gro_data["z_coord"][j]*10, 
                                _charge))

        # Write the final blank line required by Gaussian
        fout.write("\n")

    # Yield final messages
    print(f"Gaussian input file '{output}' successfully created.\n")
    if cutoff != 0.0:
        print(f"Atoms {to_remove} were removed and the charges were redistributed over the point charges within a {cutoff} AA radius sphere.\n")
    else:
        print(f"Atoms {to_remove} were removed and the charges were redistributed over the respective {n_neighbors} closest point charges.\n")


def get_distant_charge(file: str) -> float | str:
    """Get the most distant point charge from the QM region.

    Args:
        file (str): name of the Gaussian file.
    Returns:
        larg_dist (float): distance of the most distant charge.
        larg_line (str): line of the most distant charge.
    """

    # Initialize the variables
    qm_atoms = []
    larg_dist = 0
    larg_line = ''

    with open(file, "r") as f:
        line = f.readline()

        # Read lines until the second blank line is found
        i = 0
        while i < 3:
            line = f.readline()
            if line.strip() == "":
                i += 2

        # Skip the line with charge and multiplicity
        line = f.readline()
        line = f.readline()

        # Get the QM atoms coordinates       
        while line.strip() != "":
            words = line.split()

            qm_atoms.append(
                (float(words[1]),
                 float(words[2]),
                 float(words[3]))
                )
            
            line = f.readline()

        # Skip the next blank line
        line = f.readline()

        # Loop over all point charges
        while line.strip() != "":
            words = line.split()
            
            x = float(words[0])
            y = float(words[1])
            z = float(words[2])

            # Get the smallest distance to the QM atoms
            dist = 1_000_000
            for atom in qm_atoms:
                # Get the coordinates
                x_ = atom[0]
                y_ = atom[1]
                z_ = atom[2]
                
                # Compute the distance
                dist_ = ((x-x_)**2+(y-y_)**2+(z-z_)**2)**(1/2)
                if dist_ < dist:
                    dist = dist_
            
            if dist > larg_dist:
                larg_dist = dist
                larg_line = line

            line = f.readline()

    return larg_dist, larg_line


def check_total_charge(
    file: str, 
    configs_dir: str, 
    test: bool
) -> None:
    """Add the spurious charge to the most distant point charge.  

    Args:
        file (str): name of the Gaussian file.
        configs_dir (str): directory with files
        test (bool): avoid searching for charges during visualization tests.
    Returns:
        A "neutral_{file}" file.
    """

    with open(file, "r") as f:
        line = f.readline()

        if test:
            print("No partial charges when --test is set.\n")
        else:
            # Read lines until the third blank line is found
            i = 0
            while i < 3:
                line = f.readline()
                if line.strip() == "":
                    i += 1

        # Initialize the total charge
        total_charge = 0.0

        # Loop over the point charges
        while line:
            line = f.readline()
            
            # Get the total charge
            if len(line.split()) == 4:
                total_charge += float(line.split()[3])

    # Neutralize the total charge if it isn't zero
    if not test:
        if round(total_charge, 5) == 0:
            print(f"Total charge is 0.00000. No need for fixing numerical fluctuations.\n")
        else:
            # Print a warning for big charges
            if total_charge > 0.1:
                print("WARNING: Total charge is bigger than usual numerical fluctuations. " \
                      "Make sure that the charge is consistent with the expected total charge.\n")

            # Find the most distant charge (in average) to the molecule
            larg_dist, larg_line = get_distant_charge(file)

            # Write the neutralized file
            with open(file, "r") as f:
                lines = f.readlines()

                # Update the charge
                for i in range(len(lines)):
                    if lines[i] == larg_line:
                        words = lines[i].split()
                        words[-1] = str(float(words[-1]) - float(total_charge))
                        lines[i] = '\t'.join(words) + "\n"

            if configs_dir:
                # Get the name of the file
                _file = file.split("/")[-1]
                with open(f"input_files/neutral_{_file[6:]}", "w") as fout:
                    for line in lines:
                        fout.write(line)
            else:
                with open(f"neutral_{file}", "w") as fout:
                    for line in lines:
                        fout.write(line)

            print(f"Total charge is {total_charge:.5f} due to numerical fluctuations.\n"
                  f"Charge {-total_charge:.5f} was added to the most " \
                  f"distant partial charge in file 'neutral_*'.\n" \
                  f"The most distant partial charge is located at {tuple(larg_line.split()[:3])} AA " \
                  f"and the smallest distance of it to any QM atom is {larg_dist:.3f} AA.\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract configurations from Gromos87 (.gro) files.")
    
    parser.add_argument("itpfile", help="topology file(s) (.itp) from GROMACS.", nargs="+", type=str, default=[])
    parser.add_argument("--grofile", "-gro", help="reference .gro configuration file.", type=str)
    parser.add_argument("--configs_dir", "-dir", help="path to the directory with .gro configurations.", type=str, default='')
    parser.add_argument("--atomnums", "-an", help="list of atom numbers treated with QM (e.g. 1-10 22 82-93).", nargs="+", type=parse_range, default=[])
    parser.add_argument("--residues", "-res", help="list of residues treated with QM.", nargs="+", type=str, default=[])
    parser.add_argument("--resnums", "-rn", help="number of the residue(s) to be treated with QM. Default is 1.", nargs="+", type=int, default=[1])
    parser.add_argument("--link_scale_factor", "-sf", help="link atom distance scale factor. Default is 0.71.", type=float, default=0.71)
    parser.add_argument("--n_neighbors", "-nn", help="number of closest neighbors to redistribute the charge. Default is 3.", type=int, default=3)
    parser.add_argument("--cutoff", "-cut", help="cutoff radius (in AA) to redistribute the charge.", type=float, default=0.0)
    parser.add_argument("--keywords", "-k", help="Gaussian keywords for the calculation. Default is \"B3LYP/6-31G(d,p) Charge\".", nargs="+", type=str, default=["B3LYP/6-31G(d,p) Charge"])
    parser.add_argument("--charge", "-c", help="total charge of the system. Default is 0.", type=int, default=0)
    parser.add_argument("--spin_multiplicity", "-ms", help="spin multiplicity of the system. Default is 1.", type=int, default=1)
    parser.add_argument("--output", "-o", help="name of the output file. Default is \"calc.com\".", type=str, default="calc.com")
    parser.add_argument("--test", help="If True, generates a .xyz file with partial charges as bismuth atoms for visualization.", action="store_true", default=False)

    args = parser.parse_args()

    if args.grofile and args.configs_dir:
        print("Choose between a single '-gro' file or the directory 'configs_dir' with more than one file.")
        sys.exit()
    elif not args.grofile and not args.configs_dir:
        print("Provide a '-gro' file or a 'configs_dir' directory with files.")
        sys.exit()

    if args.configs_dir:
        if os.path.isdir(args.configs_dir):
            # Creates a directory for the input files
            os.makedirs("input_files", exist_ok=True)

            # Loop over all files in the configuration's directory
            for _, _, files in os.walk(args.configs_dir):
                for file in files:
                    # Track the configuration
                    _text = f"Working on configuration '{file}'..."
                    print("#" * len(_text))
                    print(f"{_text}")
                    print("#" * len(_text), "\n")

                    # Get the .gro atom numbers of QM atoms
                    qm_atoms = get_QM_atoms(os.path.join(args.configs_dir, file), args.atomnums, args.residues, args.resnums)

                    # Flatten the qm_atoms in case of nested lists
                    if any(isinstance(i, list) for i in qm_atoms):
                        qm_atoms = [atom for sublist in qm_atoms for atom in sublist]

                    # Get the .gro atom numbers of atoms to be removed and the link atoms coordinates
                    to_remove_num, link_coords = h_link_saturation(
                        os.path.join(args.configs_dir, file), 
                        args.itpfile, 
                        qm_atoms,
                        args.link_scale_factor
                    )

                    # Set a variable for the input file's name
                    _name = f"input_{file.split('.')[0]}.com"

                    write_gaussian_input(
                        os.path.join(args.configs_dir, file), 
                        args.itpfile,
                        qm_atoms,
                        to_remove_num,
                        link_coords,
                        args.n_neighbors,
                        args.cutoff,
                        args.keywords, 
                        args.charge, 
                        args.spin_multiplicity, 
                        os.path.join('input_files', _name),
                        args.test,
                    )

                    check_total_charge(
                        os.path.join('input_files', _name), 
                        args.configs_dir, 
                        args.test
                    )

        else:
            print(f"Could not access '{args.configs_dir}' directory.")

    else:
        # Get the .gro atom numbers of QM atoms
        qm_atoms = get_QM_atoms(args.grofile, args.atomnums, args.residues, args.resnums)

        # Flatten the qm_atoms in case of nested lists
        if any(isinstance(i, list) for i in qm_atoms):
            qm_atoms = [atom for sublist in qm_atoms for atom in sublist]

        # Get the .gro atom numbers of atoms to be removed and the link atoms coordinates
        to_remove_num, link_coords = h_link_saturation(
            args.grofile, 
            args.itpfile, 
            qm_atoms,
            args.link_scale_factor
        )

        write_gaussian_input(
            args.grofile, 
            args.itpfile,
            qm_atoms,
            to_remove_num,
            link_coords,
            args.n_neighbors,
            args.cutoff,
            args.keywords, 
            args.charge, 
            args.spin_multiplicity, 
            args.output,
            args.test,
        )

        check_total_charge(args.output, args.configs_dir, args.test)

        # Change the file format to .xyz
        if args.test:
            os.rename(args.output, f"{args.output[:-4]}.xyz")


if __name__ == '__main__':
    main()