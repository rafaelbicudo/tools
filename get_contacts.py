#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import contacts

import warnings
warnings.filterwarnings('ignore')


def get_contacts(u, group_a, group_b, radius=4.5, same_group=False):
    """ 
    Adapted from https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/contacts_within_cutoff.html
    """

    timeseries = []

    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)

        # Avoid couting contacts of atoms with itself
        if same_group:
            np.fill_diagonal(dist, radius+1)

        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()

        timeseries.append([ts.time, n_contacts])

    return np.array(timeseries)


def main() -> None:
    
    parser = argparse.ArgumentParser(description="Compute the number of contacts for molecules with backbone and side chains structures.")
    parser.add_argument("topfile", help="GROMACS topology/structure file path (e.g., .tpr or .gro).", type=str)
    parser.add_argument("trajfile", help="GROMACS trajectory file path (e.g. .xtc or .trr).", type=str)
    parser.add_argument("-g1", "--group1", help="Atom selection language of Group 1.", type=str, required=True)
    parser.add_argument("-g2", "--group2", help="Atom selection language of Group 2.", type=str, required=True)
    parser.add_argument("-r", "--radius", help="Cutoff radius (in Angstrom) for contacts. Default is 4.5.", type=float, default=4.5)
    parser.add_argument("-o", "--output", help="Output file name. Default is 'contacts.csv'.", type=str, default="contacts.csv")

    args = parser.parse_args()

    # Create the universe
    u = mda.Universe(args.topfile, args.trajfile)

    # Create the AtomGroups
    g1 = u.select_atoms(args.group1)
    g2 = u.select_atoms(args.group2)

    # Check for identical groups
    if g1 == g2:
        same_group=True

        print("Warning: Group 1 and Group 2 are identical. Computing contacts within the same group.")
    else:
        same_group=False

    # Compute the contacts
    conts = get_contacts(u, g1, g2, radius=args.radius, same_group=same_group)

    # Create a pandas.DataFrame and export it
    df = pd.DataFrame(conts, columns=['Time (ps)', 'Number of contacts'])
    df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
