#!/usr/bin/env python3

import os
import argparse
import sys


def get_atoms(
    file: str,
) -> int:
    """Get the number of atoms.

    Args:
        file (str): Gaussian optimization output file.

    Returns:
        n_atoms (int): number of atoms.
    """

    with open(file, "r") as f:
        line = f.readline()

        while line:
            if "NAtoms=" in line:
                n_atoms = int(line.split()[1])
                return n_atoms

            line = f.readline()


def make_traj(
    file: str,
    output: str
) -> None:
    """"Get the transient configurations and make a .xyz trajectory file with them.

    Args:
        file (str): Gaussian optimization output file.
        output (str): .xyz trajectory file name.
    """

    # Get the number of atoms
    n_atoms = get_atoms(file)

    # Open the trajectory file
    fout = open(output, "w")

    # Create a dictionary to map atom symbol and number
    Z_to_symbol = {
        1: 'H', 
        6: 'C',
        7: 'N',
        8: 'O',
        9: 'F', 
        15: 'P',
        16: 'S',
        17: 'Cl',
        35: 'Br',
        53: 'I'
    }

    with open(file, "r") as f:
        line = f.readline()

        # Initialize an index for tracking the number of optimization steps
        index = 0

        while line:
            if "Standard orientation" in line:
                index += 1

                # Write the .xyz file header
                fout.write(f'{n_atoms}\n')
                fout.write(f'Optimization step: {index}\n')

                # Ignore the Gaussian header
                for i in range(5):
                    line = f.readline()
                
                # Parse the geometry and write it in the trajectory file
                for i in range(n_atoms):
                    words = line.split()
                    
                    Z = int(words[1])
                    symbol = Z_to_symbol.get(Z)
                    x = float(words[3])
                    y = float(words[4])
                    z = float(words[5])

                    fout.write(f'\t{symbol}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n')

                    line = f.readline()

            line = f.readline()


def main() -> None:
    parser = argparse.ArgumentParser(description="Get transient configurations from geometry optimization using Gaussian.")
    
    parser.add_argument("file", help="Gaussian optimization output file.", type=str)
    parser.add_argument("-o", "--output", help="Name of the file with transient configurations. Default is opt_configs.xyz.", default="opt_configs.xyz", type=str)

    args = parser.parse_args()

    make_traj(args.file, args.output)


if __name__ == '__main__':
    main()