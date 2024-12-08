#!/usr/bin/env python3

"""Determine the glass temperature transition from GROMACS output files.

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
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


@dataclass
class Options:
    """Class to store the simulation options."""

    # Class attributes
    dirs : str = None
    vol_index: int = None
    save_df: bool = None
    vol_df: str = None
    pref: str = None
    low_temps: list = None
    high_temps: list = None
    temp_int: float = None
    t_min: float = None
    t_max: float = None
    t0: int = None
    tf: int = None

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
        "dirs": args.directories,
        "vol_index": args.volume_index,
        "save_df": args.save_data,
        "vol_df": args.volume_data,
        "pref": args.directory_preffix + '_',
        "low_temps": args.low_temperatures,
        "high_temps": args.high_temperatures,
        "temp_int": args.temperature_interval,
        "t_min": args.minimum_temperature,
        "t_max": args.maximum_temperature,
        "t0": args.initial_time,
        "tf": args.final_time
    }

    return Options(**options_dict)


def configure_run(args_in: List[str]) -> Options:
    """Configure the run based on command line arguments.

    Args:
        args_in (List[str]): The command line arguments.

    Returns:
        argparse.Namespace: The parsed command line arguments.
    """

    parser = argparse.ArgumentParser(description="Determine the glass temperature transition from GROMACS output files.")
    
    parser.add_argument(
        "-dirs",
        "--directories",
        type = str,
        nargs = "+",
        help = "the directories for the independent sets of simulation. \
                Used for more than 1 set of simulations. Default is [].",
        default = []
    )
    parser.add_argument(
        "-i",
        "--volume_index",
        type = int,
        help = "the index for volume when running 'gmx energy'. Default is 22.",
        default = 22
    )
    parser.add_argument(
        "-s",
        "--save_data",
        action = "store_true",
        help = "If True, saves the volume data in the .csv format. Default is True.",
        default = True
    )
    parser.add_argument(
        "-vd",
        "--volume_data",
        type = str,
        help = "the volume data .csv file name (e.g., vol_data.csv). Default is None.",
        default = None
    )
    parser.add_argument(
        "-pref",
        "--directory_preffix",
        type = str,
        help = "the preffix of the temperature directories. Default is 'temp'.",
        default = 'temp'
    )
    parser.add_argument(
        "-lT",
        "--low_temperatures",
        type = str,
        help = "list with low temperature values to be included in the fitting. Default is 20-70.",
        default = '20-70'
    )
    parser.add_argument(
        "-hT",
        "--high_temperatures",
        type = str,
        help = "interval of high temperature values to be included in the fitting. Default is 230-280.",
        default = '230-280'
    )
    parser.add_argument(
        "-dT", 
        "--temperature_interval", 
        type = float, 
        help = "the temperature interval between simulations (in Celsius). Default is 5.", 
        default = 5.0
    )
    parser.add_argument(
        "-T_min",
        "--minimum_temperature",
        type = float,
        help = "the minimum temperature for plotting (in Celsius). Default is 0.",
        default = 0.0
    )
    parser.add_argument(
        "-T_max",
        "--maximum_temperature",
        type = float,
        help = "the maximum temperature for plotting (in Celsius). Default is 300.",
        default = 300.0
    )
    parser.add_argument(
        "-t0",
        "--initial_time",
        type = int,
        help = "the initial time (in ps) to compute the average volume. Default is 500.",
        default = 500
    )
    parser.add_argument(
        "-tf",
        "--final_time",
        type = int,
        help = "the final time (in ps) to compute the average volume. Default is 1000.",
        default = 1000
    )

    args = parser.parse_args(args_in)

    return parse_args(args)


@dataclass
class SimulationManager:
    """Class to control the workflow."""

    options: Options

    def run_gmx_energy(self, dir, temp_dir, file):
        """Run gmx energy to get the volume data."""

        try:
            # Run the gmx energy command
            result = sp.run(
                f"gmx energy -f {os.path.join(dir, temp_dir, file)} -b {self.options.t0} -e {self.options.tf}",
                input = str(self.options.vol_index), 
                text = True,            # Ensure input/output is treated as text
                capture_output = True,  # Capture the output of the command
                check = True,           # Raise an error if the command fails
                shell = True,
            )

            # Remove the output file
            sp.run("rm energy.xvg", shell=True)

            # Find the average value in the output
            for line in result.stdout.splitlines():
                if "Volume" in line:
                    # Split the summary line
                    words = line.split()
                    
                    # Create the a data dict
                    dict_ = {
                        'dir': [temp_dir],
                        'path': [os.path.join(dir, temp_dir)],
                        'avg_volume': [float(words[1])],
                        'error': [float(words[2])],
                        'rmsd': [float(words[3])]
                    }
                
                    return dict_

        except sp.CalledProcessError as e:
            print("Error occurred while running gmx energy:", e)

    def get_volumes(self):
        """Determine the volumes."""

        # Initialize a pandas DataFrame
        df = pd.DataFrame({
            'dir': [],
            'path': [],
            'avg_volume': [],       # in nm^3
            'error': [],            # in nm^3
            'rmsd': [],             # in nm^3
            'temp': []              # in Celsius
        })

        # Loop over all independent set directories
        for dir in self.options.dirs:
            print(f"Reading data from directory: {dir}")
            # Loop over the files
            for i in os.listdir(dir):
                # Check only for the temperature directories
                if os.path.isdir(os.path.join(dir, i)) and i.startswith(self.options.pref):
                    # Loop over the temperature directories
                    for j in os.listdir(os.path.join(dir, i)):
                        # Find the energy file
                        if j.endswith('.edr'):
                            # Get the data dictionary
                            vol_dict = self.run_gmx_energy(dir, i, j)

                            # Add the temperature to the dictionary
                            vol_dict['temp'] = float(i.split('_')[-1])

                            # Add the whole data to the DataFrame
                            df = pd.concat([df, pd.DataFrame(vol_dict)], ignore_index=True)

        # Export the data as a .csv file
        if self.options.save_df:
            df.to_csv('vol_data.csv', index=False)
            print("Data exported in the 'vol_data.csv' file.\n")

        return df

    def fit_line(self, interval, df):
        """Fit a straight line to the points."""

        # Define the linear model
        def linear_model(x, a, b):
            return a * x + b

        # Initialize the arrays
        vols = []
        vol_errs = []
        temps = []

        # Loop over the temperatures
        for t in np.arange(int(interval.split('-')[0]),
                           int(interval.split('-')[1])+1,
                           self.options.temp_int):
            # Loop over the directories in the DataFrame
            for t_ in df['temp']:
                # Find the matching row with temperature
                if t == t_:
                    # Get the row data
                    row = df[df['temp'] == t_]

                    # Append the data for fitting
                    temps.append(t)
                    vols.append(row['avg_volume'].values[0])
                    vol_errs.append(row['error'].values[0])
                    
        # Fit a linear model
        popt, pcov = curve_fit(
            linear_model, 
            temps, 
            vols, 
            sigma = vol_errs, 
            absolute_sigma = True
        )

        # Get fit parameters and uncertainties
        ang_coef, lin_coef = popt
        ang_err, lin_err = np.sqrt(np.diag(pcov))

        return ang_coef, lin_coef, ang_err, lin_err


    def run(self):
        """Run the code."""

        # Get the volume data
        if self.options.vol_df:
            df = pd.read_csv(self.options.vol_df)
        else:
            df = self.get_volumes()

        # Group by temperature, drop 'path' and 'dir', and calculate the mean for other columns
        avg_df = df.drop(columns=['path', 'dir']).groupby('temp', as_index=False).mean()

        # Fit straight lines for low and high temperatures
        ang_low, lin_low, ang_low_err, lin_low_err = self.fit_line(self.options.low_temps, avg_df)
        ang_high, lin_high, ang_high_err, lin_high_err = self.fit_line(self.options.high_temps, avg_df)

        # Ensure lines are not parallel
        if ang_high == ang_low:
            raise ValueError("Lines are parallel and do not intersect. \
                              Try other combination of points")

        # Calculate x intersection
        glass_temp = (lin_high - lin_low) / (ang_low - ang_high)

        # Propagate error in x
        sigma_temp = np.sqrt(
            (lin_low_err / (ang_low - ang_high))**2 +
            (lin_high_err / (ang_low - ang_high))**2 +
            ((lin_high - lin_low) * ang_low_err / (ang_low - ang_high)**2)**2 +
            ((lin_high - lin_low) * ang_high_err / (ang_low - ang_high)**2)**2
        )

        print(f"Glass transition temperature: {glass_temp:.4f} ± {sigma_temp:.4f}\n")

        # Plot the curves
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))

        x = np.linspace(self.options.t_min, self.options.t_max, 1500)

        # Plot all points with error bars
        ax.errorbar(avg_df['temp'], avg_df['avg_volume'], yerr=avg_df['error'], fmt='o', capsize=5, color='black', alpha=0.5)

        # Get the y_values
        y_low = ang_low * x + lin_low
        y_high = ang_high * x + lin_high

        # Plot the fitted lines
        ax.plot(x, y_low, color='blue', zorder=1, label='Low temperature fit')#label = rf'$y={ang_low:.2f}x+{lin_low:.2f}$')
        ax.plot(x, y_high, color='red', zorder=2, label='High temperature fit')#label = rf'$y={ang_high:.2f}x+{lin_high:.2f}$')

        # Fill the areas of points used to fit the curves
        ax.axvspan(float(self.options.low_temps.split('-')[0]), float(self.options.low_temps.split('-')[1]), color='blue', alpha=0.3)
        ax.axvspan(float(self.options.high_temps.split('-')[0]), float(self.options.high_temps.split('-')[1]), color='red', alpha=0.3)

        # Highlight the glass transition point
        ax.scatter(glass_temp, ang_low * glass_temp + lin_low, marker='X', color='black', zorder=3, s=100)
        ax.axvline(glass_temp, 0, ang_low * glass_temp + lin_low, linestyle='--', color='black')

        ax.set_xlabel(r'Temperature ($\rm ^o$C)', fontsize=12)
        ax.set_ylabel(r'Volume (nm$^3$)', fontsize=12)
        ax.legend(loc='upper left', frameon=False, fontsize=12)
        ax.minorticks_on()
        ax.tick_params(labelsize=12, width=1.25)
        ax.grid(which='both', linestyle=':', linewidth=0.5, alpha=0.75)
        
        plt.tight_layout()
        plt.savefig('glass_temp_fit.pdf', format='pdf', dpi=600)
        plt.savefig('glass_temp_fit.png', format='png', dpi=600)
        plt.show()


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