#!/usr/bin/env python3

"""Extract excited state properties from Gaussian09/16 output file
   and export them as a .csv file.
"""

import argparse
import sys
import os


# Classes
class ExcitedState:
    def __init__(self, number: int, mult: str, energy: float, wavel: float, osc_str: float) -> None:
        """Initialize an ExcitedState instance.

        Args:
            number (int): _description_
            mult (str): _description_
            energy (float): _description_
            osc_str (float): _description_
        """

        self.number = number
        self.mult = mult
        self.energy = energy
        self.wavel = wavel
        self.osc_str = osc_str

    def info(self):
        """Display the information of the excited state."""
        print(f"Excited State {self.number}")
        print(f"Spin multiplicity: {self.mult}")
        print(f"Energy: {self.energy} eV")
        print(f"Wavelength: {self.wavel} nm")
        print(f"Oscillator Strength: {self.osc_str} arb. units")


# Functions
def get_excited_props(files: str) -> dict:
    """Extract the excited state properties from Gaussian output file.

    Args:
        file (str): Gaussian09/16 output file(s).

    Returns:
        excited_states (dict): Excited states data. 
    """

    # Initialize the variables
    exc_data = {}
    n_es = 1

    # Open and read the files
    if isinstance(files, list):
        for file in files:
            with open(file, 'r', encoding='utf-8') as f:
                lines = f.readlines()

                # Get excited state data
                for line in lines:
                    if "Excited State" in line:
                        words = line.split()

                        exc_data[n_es] = ExcitedState(
                            number = int(words[2][:-1]),
                            mult = words[3].split("-")[0],
                            energy = float(words[4]),
                            wavel = float(words[6]),
                            osc_str = float(words[8][2:])
                        )

                        n_es += 1
    else:
        with open(files, 'r', encoding='utf-8') as f:
            lines = f.readlines()

            # Get excited state data
            for line in lines:
                if "Excited State" in line:
                    words = line.split()

                    exc_data[n_es] = ExcitedState(
                        number = int(words[2][:-1]),
                        mult = words[3].split("-")[0],
                        energy = float(words[4]),
                        wavel = float(words[6]),
                        osc_str = float(words[8][2:])
                    )

                    n_es += 1
    
    return exc_data


def export_csv(logfile: str, csvfile: str) -> None:
    """Export the excited state data in the .csv format.

    Args:
        logfile (str): Gaussian09/16 output file.
        csvfile (str): Name of the .csv file.
    """

    # Get the excited states properties
    exc_data = get_excited_props(logfile)

    with open(f"{csvfile}", "w") as fout:
        # Write the header
        fout.write("number;mult;energy(eV);wavelength(nm);oscillator strength")

        # Loop over all states
        for state in exc_data:
            fout.write(f"{exc_data[state].number};"
                       f"{exc_data[state].mult};"
                       f"{exc_data[state].energy};"
                       f"{exc_data[state].wavel};"
                       f"{exc_data[state].osc_str}\n")
            

def main() -> None:
    parser = argparse.ArgumentParser(description="Extract configurations from Gromos87 (.gro) files.")
    
    parser.add_argument("-dir", help="Directory with Gaussian09/16 output files.")
    parser.add_argument("-log", "--logfile", help="Gaussian09/16 output file.", nargs="+", type=str)
    parser.add_argument( "-csv", "--csvfile", help=".csv data file name", type=str, default="exc_data.csv")

    args = parser.parse_args()

    if not args.dir and not args.logfile:
        print("Provide a single output file (-log) or a directory with more than one (-dir).\n")
        sys.exit()
    elif args.dir and args.logfile:
        print("Provide either -log or -dir.\n")
        sys.exit()

    if args.dir:
        for root, _, files in os.walk(args.dir):
            for file in files:
                if file.endswith(".log") or file.endswith(".out"):
                    export_csv(os.path.join(root, file), f"{file.split('.')[0]}-{args.csvfile}")
    else:
        export_csv(args.logfile, args.csvfile)


if __name__ == '__main__':
    main()
