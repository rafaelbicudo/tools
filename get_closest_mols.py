#!/bin/usr/env python3

"""Script to generate .xyz trajectory files with the closest 
   molecules to a given reference from GROMACS output.

   AUTHOR: Rafael Bicudo Ribeiro (@rafaelbicudo)
   DATE: 01/2024
"""

import MDAnalysis as mda
from MDAnalysis import transformations
import numpy as np
import argparse


def create_xyz_traj(
    universe,
    ref,
    molecules, 
    n_mol: int = 1,
    file_name = 'traj.xyz',
    n_frames: int = 100,
    n_init: int = 1,
    n_end: int = -1,
    cm_to_origin = False
):
    """
    Generate .xyz trajectory file with the closest 'n_mol' molecules from reference residue 'ref'.
    
    PARAMETERS:
    universe [type: MDAnalysis.core.universe.Universe] - MDAnalysis universe with all data.
    ref [type: MDAnalysis.core.groups.AtomGroup] - MDAnalysis AtomGroup with the reference molecule/residue.
    molecules [type: MDAnalysis.core.groups.AtomGroup] - MDAnalysis AtomGroup with all molecules/residues
                                                         to search for the closest 'n_mol'.
    n_mol [type: int] - Number of closest residues to the reference. 
    file_name [type: str] - Name of the output trajectory file.
    n_frames [type: int] - Number of frames considered.
    n_init [type: int] - First frame.
    n_end [type: int] - Last frame.
    cm_to_origin [type: bool] - If True, translate the target's residue center of mass to origin.
    
    OUTPUT:
    'file_name' file.
    """

    with open("{}".format(file_name), "w") as fout:

        # Get the step to get n_frames
        if n_end == -1: # For some reason this is not printing the .xyz file afterwards
            dt = round((len(universe.trajectory)-n_init)/n_frames)
        else:
            dt = round((n_end-n_init)/n_frames)

        # Get the final frame
        if n_end == -1:
            final_frame = n_end
        else:
            final_frame = n_end + 1

        # Loop over all frames
        for ts in universe.trajectory[n_init:final_frame:dt]:

            # Calculate the center of mass for each molecule (residue) in the molecule group
            molecule_centers = np.array([res.atoms.center_of_mass() for res in molecules.residues])

            # Calculate the distances between the center of mass of the reference group and each molecule
            distances = np.linalg.norm(molecule_centers - ref.center_of_mass(), axis=1)

            # Find the indices of the 'n_mol' closest molecules
            closest_indices = np.argsort(distances)[:n_mol]

            # Get the residues of the 'n_mol' closest molecules
            closest_molecules = [molecules.residues[i].atoms for i in closest_indices]
            
            # Shift the target's residue center of mass to origin
            if cm_to_origin:
                
                # Compute the center of mass of ref
                cm = ref.center_of_mass()
                
                # Translate all the system
                universe.atoms.positions -= cm
            
            # Get the number of atoms
            n_atoms = len(ref)
            
            for molecule in closest_molecules:
                n_atoms += len(molecule)
            
            # Write the .xyz trajectory file header
            fout.write("{}\n".format(n_atoms))
            fout.write("Frame {}\n".format(ts.frame))

            # Write the reference molecule coordinates
            for atom in ref:
                fout.write("{:>8}{:>12.5f}{:>11.5f}{:>11.5f}\n".format(atom.element, atom.position[0], 
                                                                       atom.position[1], atom.position[2]))

            # Write the coordinates of the 'n_mol' closest molecules
            for molecule in closest_molecules:
                for atom in molecule:
                    fout.write("{:>8}{:>12.5f}{:>11.5f}{:>11.5f}\n".format(atom.element, atom.position[0], 
                                                                           atom.position[1], atom.position[2]))


def main() -> None:
        
    parser = argparse.ArgumentParser(
        description = "Get closest molecules to a given reference from GROMACS output."
    )
    
    parser.add_argument(
        "-s",
        "--topology_file",
        help = "GROMACS structure+mass(db) file (.tpr, .gro, .g96, .pdb, .brk or .ent). Default is 'topol.tpr'.",
        default = 'topol.tpr',
        type = str
    )
    parser.add_argument(
        "-f",
        "--trajectory_file",
        help = "GROMACS trajectory file (.xtc, .trr, .cpt, .gro, .g96, .pdb or .tng). Default is 'traj.xtc'.",
        default = 'traj.xtc',
        type = str    
    )
    parser.add_argument(
        "-ref",
        "--reference", 
        help = "Atom selection language (same used in VMD) for the reference group. Default is 'resname UNK'.",
        default = 'resname UNK',
        type = str
    )
    parser.add_argument(
        "-cref",
        "--closest_mols_reference", 
        help = "Atom selection language (same used in VMD) for the closest atoms. Default is 'not resname UNK'.",
        default = 'not resname UNK',
        type = str
    )
    parser.add_argument(
        "-n",
        "--n_closest", 
        help = "Number of closest molecules to be included. Default is 1.",
        default = 1,
        type = int
    )
    parser.add_argument(
        "-init",
        "--initial_frame",
        help = "Initial frame. Default is 1.",
        default = 1,
        type = int
    )
    parser.add_argument(
        "-end",
        "--final_frame",
        help = "Final frame. Default is -1 (last frame).",
        default = -1,
        type = int
    )
    parser.add_argument(
        "-n_frames",
        "--number_of_frames",
        help = "Number of frames to be considered. Default is 100.",
        default = 100,
        type = int
    )
    parser.add_argument(
        "-cm_to_origin",
        "--center_of_mass_to_origin",
        help = "If True, translates the center of mass of the reference molecule to the origin. Default is False.",
        default = False,
        action = "store_true"
    )
    parser.add_argument(
        "-o",
        "--output_file",
        help = "Output file name. Default is 'output.xyz'.",
        default = 'output.xyz',
        type = str
    )

    args = parser.parse_args()

    u = mda.Universe(
        args.topology_file, 
        args.trajectory_file
    )

    ref_atoms = u.select_atoms(args.reference)
    closest_mols = u.select_atoms(args.closest_mols_reference)

    # Centralize the reference atoms and avoid pbc problems
    # https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/
    workflow = (transformations.unwrap(u.atoms),
                    transformations.center_in_box(ref_atoms, center='mass'),
                    transformations.wrap(u.atoms, compound='fragments'))
    u.trajectory.add_transformations(*workflow)

    create_xyz_traj(
        universe = u, 
        ref = ref_atoms, 
        molecules = closest_mols, 
        n_mol = args.n_closest, 
        file_name = args.output_file, 
        n_frames = args.number_of_frames, 
        n_init = args.initial_frame, 
        n_end = args.final_frame, 
        cm_to_origin = args.center_of_mass_to_origin
    )

    print("\nTrajectory file '{}' has been created.".format(args.output_file))

if __name__ == "__main__":
    main()





