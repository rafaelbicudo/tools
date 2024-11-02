#!/usr/bin/env python3

import os
import argparse
import sys
import MDAnalysis as mda
import numpy as np
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation as R


def stack(
    file: str,
    type: str,
    dist: float,
    output: str,
    shear: float,
    rot_angle: float
):
    """Create an interface.

    Args:
        file (str): file with configuration.
        type (str): type of aggregation.
        dist (float): distance between molecules.
        output (str): name of the output file.
        shear (float): distance to be sheared about (for J-type).
        rot_angle (float): rotation angle in degrees (for T- and X-type).
    """

    # Read the file
    u = mda.Universe(file)

    # Performs a PCA on the coordinates
    pca = PCA(n_components=3)
    pca.fit(u.atoms.positions)

    # The third principal component is usually normal
    # to the molecule's plane
    pca3 = pca.components_[-1]

    if type == "H":
        # Get the coordinates for the new shifted molecule
        new_coords = u.atoms.positions + pca3 * dist

    elif type == "J":
        # Get the principal component axis
        pca1 = pca.components_[0]

        # Get the coordinates for the new shifted and sheared molecule
        new_coords = u.atoms.positions + pca3 * dist + pca1 * shear

    elif type == "X":
        # Compute the rotation angle in radians
        angle = np.pi * rot_angle / 180

        # Get the rotated coordinates
        rotation = R.from_rotvec(angle * pca3)
        new_coords = rotation.apply(u.atoms.positions)

        # Get the coordinates for the new shifted molecule
        new_coords = new_coords + pca3 * dist

    elif type == "T":
        # Get the principal component axis
        pca1 = pca.components_[0]

        # Compute the rotation angle in radians
        angle = np.pi * rot_angle / 180

        # Get the rotated coordinates
        rotation = R.from_rotvec(angle * pca1)
        new_coords = rotation.apply(u.atoms.positions)

        # Get the coordinates for the new shifted molecule
        new_coords = new_coords + pca3 * dist

    else:
        print(f"{type}-type aggregation is not available. Only 'H', 'J', 'X' or 'T' types are available.")
        sys.exit()

    # Create an universe with duplicated atoms
    new_u = mda.Merge(u.atoms, u.atoms)

    # Assign shifted coordinates to the new atoms
    new_u.atoms[u.atoms.n_atoms:].positions = new_coords

    # Export the new configuration
    new_u.atoms.write(output)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create interface between molecules.")
    
    parser.add_argument("file", help="File with configuration parsed by MDAnalysis.", type=str)
    parser.add_argument("-t", "--type", help="Type of aggregation, i.e., H-, X-, T- or J-type. Default is H.", default="H", type=str)
    parser.add_argument("-d", "--distance", help="Distance between molecules. Default is 3 Angstrom.", default=3.0, type=float)
    parser.add_argument("-o", "--output", help="Name and format of the output file. Default is stacked_file.pdb.", default="stacked_file.pdb", type=str)
    parser.add_argument("-s", "--shear", help="Shear distance to shift the molecule for J-type aggregation. Default is 1.5 Angstrom.", default=1.5, type=float)
    parser.add_argument("-rot", "--rot_angle", help="Rotation angle for the X- or T-type aggregation. Default is 90 degrees.", default=90, type=float)
    # parser.add_argument("-n", "--number", help="Number of molecules to be stacked. Default is 2.", default=2, type=int)

    args = parser.parse_args()

    stack(
        args.file, 
        args.type, 
        args.distance, 
        args.output, 
        args.shear, 
        args.rot_angle
        # args.number
    )


if __name__ == '__main__':
    main()


