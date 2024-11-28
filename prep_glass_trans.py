#!/usr/bin/env python3

"""Generate the input files for running the batch of simulations for
   determining the glass transition temperature.

   AUTHOR: Rafael Bicudo Ribeiro (@rafaelbicudo)
   DATE: 11/2024
"""

from dataclasses import dataclass
from typing import List
import time
import os
import sys
import argparse
import numpy as np
import subprocess as sp


@dataclass
class Options:
    """Class to store the simulation options."""

    # Class attributes
    gro_file: str = None
    top_file: str = None
    init_temp: float = None
    final_temp: float = None
    temp_int: float = None
    sub_file: str = None

    # Functions
    def update(self, new: dict) -> None:
        """Update the instance with new values.

        Args:
            new (dict): A dictionary containing the new values to update.

        Returns:
            None
        """
        for key, value in new.items():
            if hasattr(self, key):
                setattr(self, key, value)


def parse_args(args: argparse.Namespace) -> Options:
    """Parse all the information from the argument parser, storing in the
    SimulationOptions class.

    Define file names and set the pointers to the correct functions.

    Args:
        args (argparse.Namespace): The arguments parsed from the argument 
                                   parser.

    Returns:
        Options: An instance of the SimulationOptions class with the 
                        parsed options.
    """

    options_dict = {
        "gro_file": args.gro_file,
        "top_file": args.top_file,
        "init_temp": args.initial_temperature,
        "final_temp": args.final_temperature,
        "temp_int": args.temperature_interval,
        "sub_file": args.submission_file
    }

    return Options(**options_dict)


def configure_run(args_in: List[str]) -> Options:
    """Configure the run based on command line arguments.

    Args:
        args_in (List[str]): The command line arguments.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """

    parser = argparse.ArgumentParser(description="Prepare GROMACS input files for the simulations required to determine the glass transition temperature.")
    
    parser.add_argument(
        "gro_file",
        type = str,
        help = "the .gro input file."
    )
    parser.add_argument(
        "top_file",
        type = str,
        help = "the .top topology file."
    )
    parser.add_argument(
        "-T0", 
        "--initial_temperature", 
        type = float, 
        help = "the initial temperature (in Celsius). Default is 20.", 
        default = 20.0
    )
    parser.add_argument(
        "-Tf", 
        "--final_temperature", 
        type = float, 
        help = "the final temperature (in Celsius). Default is 280.", 
        default = 280.0
    )
    parser.add_argument(
        "-dT", 
        "--temperature_interval", 
        type = float, 
        help = "the temperature interval between simulations (in Celsius). Default is 5.", 
        default = 5.0
    )
    parser.add_argument(
        "-sub",
        "--submission_file",
        type = str,
        default = 'sub.txt',
        help = "file with the header of the submission file. Default is sub.txt."
    )

    args = parser.parse_args(args_in)

    return parse_args(args)


@dataclass
class SimulationManager:
    """Class to control the workflow."""

    options: Options
    
    def prep_sub_file(self, t):
        """Write the submission file."""
    
        # Create the submission file
        fout = open(f'run_temp{int(t)}', "w")

        # Read the header file
        with open(self.options.sub_file, "r") as f:
            lines = f.readlines()
            
            # Write the header into the submission file
            for line in lines:
                fout.write(line)

        # Move the file to the current directory
        sp.run(f"mv run_temp{int(t)} temp_{int(t)}", shell=True)

    def change_temp(self, t):
        """Change the temperature in the .mdp files and 
           move them to the current directory."""

        # Loop over the files
        for file in os.listdir(os.getcwd()):

            # Change the reference temperature in the .mdp file 
            if file.endswith('mdp'):
                with open(file, "r") as f:
                    lines = f.readlines()
                
                # Initialize the variables
                found_ref_t = False
                new_lines = []

                # Loop over the files
                for line in lines:
                    # Change the reference temperature (in Kelvin)
                    if "ref_t" in line:
                        new_lines.append(f";{line.strip()}\n")
                        new_lines.append(f"ref_t \t\t\t= {t+273}\n")
                        found_ref_t = True
                    else:
                        new_lines.append(line)

                # Warning message if temperature was not found
                if not found_ref_t:
                    print(f"Couldn't find the reference temperature in '{file}' file. \n")

                # Write the new .mdp file
                with open(os.path.join(f'temp_{int(t)}', file), "w") as fout:
                    fout.writelines(new_lines)

    def run_all(self):
        """Create a script to submit all files at once."""

        # List all directories in the current working directory
        directories = [d for d in os.listdir() if os.path.isdir(d) and d.startswith("temp_")]

        # Start building the script content
        script_content = "#!/bin/bash\n\n"
        for directory in directories:
            # Extract the number part (X) from "temp_X"
            temp_number = directory.split('_')[1]
            
            # Add the command for this directory to the script
            script_content += f"if [ -f {directory}/run_temp{temp_number} ]; then\n"
            script_content += f"  (cd {directory} && sbatch run_temp{temp_number})\n"
            script_content += "fi\n\n"

        # Write the content to `runall.sh`
        with open("runall.sh", "w") as script_file:
            script_file.write(script_content)

        # Make the script executable
        os.chmod("runall.sh", 0o755)

        print("runall.sh has been created.\n")

    def run(self):
        """Run the code."""

        # Loop over all the temperatures
        for t in np.arange(self.options.init_temp, self.options.final_temp+1, self.options.temp_int):
            
            # Create a directory for the current temperature
            sp.run(['mkdir', f'temp_{int(t)}'], shell=False)

            # Copy input files to the current directory
            sp.run(f'cp *itp {self.options.gro_file} {self.options.top_file} temp_{int(t)}', shell=True)

            # Modify and copy the .mdp files
            self.change_temp(t)

            # Create the submission file
            self.prep_sub_file(t)
            
        # Create a bash script to submit all files at once
        self.run_all()


def main(args: List[str] = None) -> None:

    # Start the time count
    global_start_time = time.monotonic()

    # Parse command-line arguments
    if args is None:
        args = sys.argv[1:]

    # Get the classes with parsed arguments
    options = configure_run(args)

    manager = SimulationManager(
        options = options
    )

    manager.run()

    # Stop the clock
    global_end_time = time.monotonic()

    # Print the time for running the code
    print(f"Total time: {global_end_time - global_start_time:.6f} sec.\n")


if __name__ == '__main__':
    main()