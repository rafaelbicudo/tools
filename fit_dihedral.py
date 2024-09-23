#!/usr/bin/env python3 

"""
Script do the following:

1) Extract QM energy from torsional scan via Gaussian09/16 and write a "qm_scan.csv" file 
2) Determine all dihedral angles that change during the scan and write a "dihedrals.csv" file 
3) Interpolate the QM data using a cubic polynomial and perform a linear regression to determine the 6 coefficients
that yielded better results.
4) Write the 6 coefficients in the topology "LR_*.itp" file.

[ ] Create a function that searchs for outliers and remove them from the linear regression.
[X] Create an option to increase the relevance of mininum points, similar to DICEtools/fit_torsional.py.
[X] Fix the linear fit problems when removing a few angles during the scan due to the atomic nuclei overlapping. 
	It seems that distances smaller than 0.5 Angstrom crash the code.
[ ] Add an option to determine the RB coefficients from the total energy instead of the equivalent "QM" torsional energy.
[ ] Add an option to read a list of angles and fit the coefficients only to these dihedrals.

Author: Rafael Bicudo Ribeiro (@rafaelbicudo) and Thiago Duarte (@thiagodsd)
DATE: DEZ/2022
"""

import argparse
import numpy as np
import os

from plot_en_angle_gaussian_scan import parse_en_log_gaussian
from plot_eff_tors import get_phi, get_potential_curve
from parse_gaussian_charges import find_natoms
from fit_torsional import shift_angle_rad

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model  import LinearRegression, Ridge, Lasso, LassoCV
from scipy.interpolate import CubicSpline, interp1d


def find_bonded_atoms(topfile: str, a1: int, a2: int, a3: int, a4: int) -> list:
	"""
	Return a list of all atoms which are candidates to change during the rotation of a1-a2-a3-a4 angle.
	
	PARAMETERS:
	topfile [type: str] - topology (.itp) file
	a1 [type: int] - first atom defining the dihedral angle
	a2 [type: int] - second atom defining the dihedral angle
	a3 [type: int] - third atom defining the dihedral angle
	a4 [type: int] - fourth atom defining the dihedral angle

	OUTPUT:
	Angles [type: list]
	"""

	dih_atoms = [a1, a2, a3, a4]
	dih_atoms = [int(i) for i in dih_atoms]

	_atoms = dih_atoms.copy()

	with open(topfile, "r") as f:
		line = f.readline()
		while "[ bonds ]" not in line:
			line = f.readline()

		line = f.readline()
		while (len(line.split()) != 0) and not (line.strip().startswith("[")):
			words = line.split()

			if (int(words[0]) in dih_atoms) and (int(words[1]) not in dih_atoms):
				_atoms.append(int(words[1]))

			elif (int(words[1]) in dih_atoms) and (int(words[0]) not in dih_atoms):
				_atoms.append(int(words[0]))
			line = f.readline()

	return _atoms


def find_dihedrals(topfile: str, a1: int, a2: int, a3: int, a4: int) -> list:
	"""
	Find which torsionals are candidates for changing during the scan.

	PARAMETERS:
	topfile [type: str] - topology (.itp) file
	a1 [type: int] - first atom defining the dihedral angle
	a2 [type: int] - second atom defining the dihedral angle
	a3 [type: int] - third atom defining the dihedral angle
	a4 [type: int] - fourth atom defining the dihedral angle

	OUTPUT: 
	Dihedral angles [type: list] (list of lists with 4 integers each)
	Labels [type: list] (list of strings)
	"""

	candidates = find_bonded_atoms(topfile, a1, a2, a3, a4)

	dih_atoms = [a1, a2, a3, a4]
	dih_atoms = [int(i) for i in dih_atoms]

	tors = []
	tors_label = []

	with open(topfile, "r") as f:
		line = f.readline()
		while "[ dihedrals ]" not in line:
			line = f.readline()

		while True:
			line = f.readline()
			if line.strip().startswith(";"):
				continue
			if (len(line.split()) == 11) and (int(line.split()[4]) != 2) and (int(line.split()[4]) != 4):
				break
			if not line:
				print("There are no proper dihedrals.")
				break

		while (len(line.split()) != 0) and not (line.strip().startswith("[")):
			words = line.split()
			_atoms = [words[0], words[1], words[2], words[3]]
			_atoms = [int(i) for i in _atoms]
			if all(elem in candidates for elem in _atoms):
				tors_label.append("%s-%s-%s-%s" % (words[0], words[1], words[2], words[3]))
				tors.append([int(words[0]), int(words[1]), int(words[2]), int(words[3])])
			line = f.readline()

	return tors, tors_label


def get_dihedrals(xyzrotationsfile: str, topfile: str, a1: int, a2: int, a3: int, a4: int, npoints: int, angles: list = []) -> pd.DataFrame:
	"""
	Read the .xyz file and find torsional angle changes.

	PARAMETERS:
	xyzrotationsfile [type: str] - configurations file
	topfile [type: str] - topology (.itp) file
	a1 [type: int] - first atom defining the dihedral angle
	a2 [type: int] - second atom defining the dihedral angle
	a3 [type: int] - third atom defining the dihedral angle
	a4 [type: int] - fourth atom defining the dihedral angle
	npoints [type: int] - number of configurations during the scan
	angles [type: np.ndarray] - list with angles to be used in the fit (default is to use all angles)

	OUTPUT:
	df_fit [type: pd.DataFrame] - dihedral angles which changed during the scan and will be used for fitting
	df [type: pd.DataFrame] - dihedral angles which changed during the scan
	"""

	dih = []
	atomsCoord = {}
	i = 0

	# Get the number of atoms in the molecule
	natoms = find_natoms(topfile)

	# Get the torsional angles and labels
	tors, tors_label = find_dihedrals(topfile, a1, a2, a3, a4)

	# Create the data frame with proper size
	df = pd.DataFrame(np.zeros((npoints, len(tors))), index=range(1, npoints+1), columns=tors_label)

	with open(xyzrotationsfile) as xyz_f:
		line = xyz_f.readline()
		words = line.split()

		# Loop over all configurations
		for conf in range(npoints):
			while "Dihedral = " not in line:
				line = xyz_f.readline()
			
			dih.append(float(line.split()[2]))

			while len(words) != 4:
				line = xyz_f.readline()
				words = line.split()

			# Parse the atomic coordinates to a dictionary (adapted from "parse_txt" function of DICEtools)
			anum = 1
			for i in range(natoms):
				atomsCoord[anum] = [float(x) for x in line.split()[1:4]]
				anum += 1
				line = xyz_f.readline()
				words = line.split()

			# Add the torsional angles to the dataframe
			for _list in tors:
				df.iat[conf, tors.index(_list)] = get_phi(atomsCoord[_list[0]], atomsCoord[_list[1]], atomsCoord[_list[2]], atomsCoord[_list[3]])
			conf += 1

		for label in tors_label:
			column = round(df[label], 1).to_numpy()
			if (column == column[0]).all():
				del df[label]

	# Add the dihedral angle to the dataframe
	dih = np.array(dih)*np.pi/180.
	df.insert(0, 'Dihedral (rad)', dih)

	# Sort the df following the dihedral angle
	df = df.sort_values(by='Dihedral (rad)', ascending=True)
	df = df.reset_index(drop=True)

	# Select only the angles indicated to be fitted
	if angles:
		# Convert the angles to a numpy array and change it to radians
		angles = np.asarray(angles)
		angles = angles * np.pi/180
		
		# Create a list with indexes to be considered
		fit = []

		# Find and remove the angles not required
		for i in range(len(df['Dihedral (rad)'])):
			for j in angles:
				if np.abs(df.loc[i, 'Dihedral (rad)'] - j) < 0.08:
					fit.append(i)

		df_fit = df.loc[df.index.isin(fit)]
	else:
		df_fit = df

	return df_fit, df


def check_overlap(topfile: str, xyzrotationsfile: str, cutoff: float, npoints: int) -> list:
	"""
	Read the nuclei positions to search for atomic overlap.

	PARAMETERS:
	topfile [type: str] - topology (.itp) file
	xyzrotationsfile [type: str] - configurations file
	cutoff [type: float] - minimum atomic distance tolerated
	npoints [type: int] - number of configurations during the scan

	OUTPUT: 
	Prints the overlapping atoms at a given dihedral configuration.
	overlap_dih [type: list[float]] - list with overlapping dihedrals
	conf_dih [type: list[float]] - list with overlapping dihedrals configurations 
	"""

	overlap_dih = []
	conf_dih = []
	atomsCoord = {}
	i = 0

	# Get the number of atoms in the molecule
	natoms = find_natoms(topfile)

	with open(xyzrotationsfile) as xyz_f:
		line = xyz_f.readline()
		words = line.split()

		# Loop over all configurations
		for conf in range(npoints):
			while len(words) != 4:
				line = xyz_f.readline()
				words = line.split()

				if "Dihedral = " in line:
					dih = words[2]

			# Parse the atomic coordinates to a dictionary (adapted from "parse_txt" function from DICEtools)
			anum = 1
			for i in range(natoms):
				atomsCoord[anum] = [float(x) for x in line.split()[1:4]]
				anum += 1
				line = xyz_f.readline()
				words = line.split()

			# Loop over atomic pairs
			j = 1
			for j in range(1, natoms):
				for i in range(j+1, natoms):
					# Compute the distance between atoms
					dist = np.linalg.norm(np.array(atomsCoord[i]) - np.array(atomsCoord[j]))

					# Print a warning if overlapping and add the dihedral angle to the 'overlap_dih' list.
					if dist <= cutoff:
						if __name__ == "__main__":
							print("The distance between atoms {} and {} from \'Dihedral = {}\' (conf = {}) is {} Angs.\n"
								"Run --remove-overlap flag to remove the overlapping configurations.\n".format(j, i, dih, conf, round(dist, 3)))
						else:
							print("Removing \'Dihedral = {}\' (conf = {}) due to atomic overlap.".format(dih, conf))
						
						if float(dih) not in overlap_dih:
							overlap_dih.append(float(dih)) 
							conf_dih.append(conf)

	return overlap_dih, conf_dih


def remove_overlap(scan: pd.DataFrame, dih_fit: pd.DataFrame, topfile: str, xyzrotationsfile: str, cutoff: float, npoints: int) -> None:
	"""
	Remove the configurations with atomic overlap.

	PARAMETERS:
	scan [type: pd.DataFrame] - data from classical and quantum scan simulations
	topfile [type: str] - topology (.itp) file
	xyzrotationsfile [type: str] - configurations file
	cutoff [type: float] - minimum atomic distance tolerated
	npoints [type: int] - number of configurations during the scan

	OUTPUT:
	Files .csv without the overlapping configurations.
	"""

	# Get the overlapping configurations
	overlap_dih, conf_dih = check_overlap(topfile, xyzrotationsfile, cutoff, npoints)

	for conf in conf_dih:
		scan = scan.drop(conf-1)
		dih_fit = dih_fit.drop(conf-1)

	scan = scan.reset_index(drop=True)
	dih_fit = dih_fit.reset_index(drop=True)

	return scan, dih_fit


def get_data(gaussianlogfile: str, txtfile: str, dfrfile: str, a1: int, a2: int, a3: int, a4: int) -> pd.DataFrame:
	"""
	Read the gaussian output file and write the .csv file with dihedrals and energies.

	PARAMETERS:
	gaussianlogfile [type: str] - output file from Gaussian09/16
	txtfile [type: str] - DICE .txt file
	dfrfile [type: str] - DICE .dfr file
	a1 [type: int] - first atom defining the dihedral angle
	a2 [type: int] - second atom defining the dihedral angle
	a3 [type: int] - third atom defining the dihedral angle
	a4 [type: int] - fourth atom defining the dihedral angle

	OUTPUT:
	scan [type: pd.DataFrame] - data from classical and quantum scan simulations
	scan_fit [type: pd.DataFrame] - data from classical and quantum scan simulations to be fitted
	"""

	# Read the dihedrals (in degrees) and energies (in kcal/mol)
	died, enqm = parse_en_log_gaussian(gaussianlogfile)

	# Check for repeated dihedrals and delete them
	for i in range(1, len(died)-1):
		if died[i] - died[i-1] < 0.001:
			del died[i]
			del enqm[i]

	# Change dihedrals to rad
	died = [shift_angle_rad(x*np.pi/180.) for x in died]

	# Extract the non-bonded potentials and classical dihedrals
	diedClass, diedEn, nben, _ = get_potential_curve(txtfile, dfrfile, a1, a2, a3, a4, died, "", False, False, False, False)

	# Consistency check for dihedral angles in classical and quantum calculations
	diff = [died[i]-diedClass[i] for i in range(len(died))]

	for angle in diff:
		if angle > 0.001:
			print("Dihedral angles do not follow the same order, please check .log and .xyz files")
			print("Quantum dihedral - Classical dihedral = ", diff)
			exit()

	# Convert lists to numpy arrays
	died = np.asarray(died)
	enqm = np.asarray(enqm)
	nben = np.asarray(nben)
	diedEn = np.asarray(diedEn)

	# Create a pd.DataFrame
	scan = pd.DataFrame(
		{
		'Dihedral (rad)': died,
		'E_qm (kcal/mol)': enqm,
		'Non-bonded (kcal/mol)': nben, 
		}
	)

	return scan


def rescale_data(scan: pd.DataFrame, angles: list = [], maxvalue=None) -> pd.DataFrame:
	"""
	Rescale the data to perform the linear regression.

	PARAMETERS:
	scan [type: pd.DataFrame] - data from classical and quantum scan simulations
	angles [type: np.ndarray] - list with angles to be used in the fit (default is to use all angles)
	maxvalue [type: float] - truncate the energy up to [-maxvalue, maxvalue]
	"""

	# Select only the angles indicated to be fitted
	if angles:
		# Convert the angles to a numpy array and change it to radians
		angles = np.asarray(angles)
		angles = angles * np.pi/180
		
		# Create a list with indexes to be considered
		fit = []

		# Find and remove the angles not required
		for i in range(len(scan['Dihedral (rad)'])):
			for j in angles:
				if np.abs(scan.loc[i, 'Dihedral (rad)'] - j) < 0.08:
					fit.append(i)

		scan_fit = scan.loc[scan.index.isin(fit)]

	else:
		scan_fit = scan

	# Create a copy to avoid warning messages
	scan_fit = scan_fit.copy()

	# Subtract the non-bonded energy from lower QM energy configuration
	enqm_0 = scan_fit.loc[:, 'E_qm (kcal/mol)'].min()
	nben_0 = scan_fit.loc[scan_fit.loc[:, 'E_qm (kcal/mol)'] == enqm_0, 'Non-bonded (kcal/mol)'].values[0]

	# Subtract the non-bonded energy from both scans
	scan_fit.loc[:, 'Non-bonded (kcal/mol)'] = scan_fit.loc[:, 'Non-bonded (kcal/mol)'] - nben_0
	scan.loc[:, 'Non-bonded (kcal/mol)'] = scan.loc[:, 'Non-bonded (kcal/mol)'] - nben_0

	# Rescale the QM energy
	scan_fit.loc[:, 'E_qm (kcal/mol)'] = scan_fit.loc[:, 'E_qm (kcal/mol)'] - enqm_0
	scan.loc[:, 'E_qm (kcal/mol)'] = scan.loc[:, 'E_qm (kcal/mol)'] - enqm_0

	# Determine the "QM" torsional energy
	scan_fit.loc[:, 'U_tors-qm (kcal/mol)'] =  scan_fit.loc[:, 'E_qm (kcal/mol)'] - scan_fit.loc[:, 'Non-bonded (kcal/mol)']
	scan.loc[:, 'U_tors-qm (kcal/mol)'] = scan.loc[:, 'E_qm (kcal/mol)'] - scan.loc[:, 'Non-bonded (kcal/mol)']

	# Rename the columns
	scan_fit = scan_fit.rename(columns={'Non-bonded (kcal/mol)': 'Non-bonded - nben0 (kcal/mol)', 
										'E_qm (kcal/mol)': 'E_qm - enqm_0 (kcal/mol)'})
	scan = scan.rename(columns={'Non-bonded (kcal/mol)': 'Non-bonded - nben0 (kcal/mol)', 
								'E_qm (kcal/mol)': 'E_qm - enqm_0 (kcal/mol)'})

	# Check for energy barrier limit request
	if maxvalue is not None:
		# enfit = np.clip(enfit, -float(maxvalue), float(maxvalue))

		# Clipping a specific column
		scan_fit['U_tors-qm (kcal/mol)'] = scan_fit['U_tors-qm (kcal/mol)'].clip(lower=-float(maxvalue), upper=float(maxvalue))
		scan['U_tors-qm (kcal/mol)'] = scan['U_tors-qm (kcal/mol)'].clip(lower=-float(maxvalue), upper=float(maxvalue))

	return scan_fit, scan


def write_itp_file(topfile: str, lr_data: dict):
	"""
	Write the linear regression coefficients in the topology file.

	PARAMETERS:
	topfile [type: str] - topology (.itp) file
	lr_data [type: dict] - nested dictionary with atoms and coeficients for each torsional angle

	OUTPUT:
	A "LR_*.itp" file.
	"""

	fout = open("LR_" + os.path.basename(topfile), "w")
	fout.write("; ########################################################################### \n")
	fout.write("; ### LINEAR REGRESSION WAS PERFORMED TO DETERMINE TORSIONAL COEFFICIENTS ### \n")
	fout.write("; ########################################################################### \n")

	with open(topfile, "r") as f:
		line = f.readline()

		while "[ dihedrals ]" not in line:
			fout.write(line)
			line = f.readline()

		while "[ pairs ]" not in line:
			fout.write(line)
			line = f.readline()

			# Change the data from kcal/mol to kJ/mol and parse it to the topology file.
			if len(line.split()) == 11:
				words = line.split()
				if words[4] == '3':
					_tors = [words[0], words[1], words[2], words[3]]
					for k in range(len(lr_data)):
						if all(elem in lr_data[k]["atoms"] for elem in _tors):
							line = line.replace(words[5], str(round(lr_data[k]["constants"][0]*4.184, 3)), 1)
							line = line.replace(words[6], str(round(lr_data[k]["constants"][1]*4.184, 3)), 1)
							line = line.replace(words[7], str(round(lr_data[k]["constants"][2]*4.184, 3)), 1)
							line = line.replace(words[8], str(round(lr_data[k]["constants"][3]*4.184, 3)), 1)
							line = line.replace(words[9], str(round(lr_data[k]["constants"][4]*4.184, 3)), 1)
							line = line.replace(words[10], str(round(lr_data[k]["constants"][5]*4.184, 3)), 1)
				else:
					print("Not a Ryckaert-Bellemans dihedral.")

		while line:
			fout.write(line)
			line = f.readline()


def linear_regression(scan: str, scan_fit: str, dih_fit: str, topfile: str, method: str, lasso_alpha: float, weight_minimums: float, fit_from_total: bool) -> None:
	"""
	Perform the linear regression, generate a plot with classical and "quantum" torsional energies
	and change the C_i coefficients in the topology file.

	PARAMETERS:
	scan [type: pd.DataFrame] - data from classical and quantum scan simulations
	scan_fit [type: pd.DataFrame] - data from classical and quantum scan simulations to be fitted
	dih [type: pd.DataFrame] - dihedral angles which changed during the scan
	dih_fit [type: pd.DataFrame] - dihedral angles which changed during the scan and will be used for fitting
	topfile [type: str] - topology (.itp) file
	method [type: str] - the method used in the linear regression (least-square, ridge or lasso)
	weight_minima [type: list] - the weights atributted to the total energy minimum points
	fit_from_total [type: bool] - fit the dihedral from total energy

	OUTPUT:
	The "linear_regression.png" and "LR_*.itp" files.
	"""

	if fit_from_total:
		print("Not implemented yet.")
		exit()

	else:
		# Cubic polynomial interpolation
		interpolarizer = CubicSpline(scan_fit.loc[:, 'Dihedral (rad)'], scan_fit.loc[:, 'U_tors-qm (kcal/mol)'])

	# Add the interpolation in the dataFrame data
	dih_fit = dih_fit.copy()
	dih_fit.loc[:, "ans"] = dih_fit.apply(lambda x: interpolarizer(x['Dihedral (rad)']), axis=1)

	# Create a new data frame with ans values
	df = dih_fit[["ans"]].copy().rename(columns={"ans":"y"})

	# Compute each Ryckaert-Bellemans torsional element
	for c in dih_fit.drop(columns=["ans", "Dihedral (rad)"]).columns:
		df[f"{c}_C0"] = np.cos(dih_fit[c])**0
		df[f"{c}_C1"] = np.cos(dih_fit[c])**1
		df[f"{c}_C2"] = np.cos(dih_fit[c])**2
		df[f"{c}_C3"] = np.cos(dih_fit[c])**3
		df[f"{c}_C4"] = np.cos(dih_fit[c])**4
		df[f"{c}_C5"] = np.cos(dih_fit[c])**5

	X = df.drop(columns=["y"])
	y = df["y"]

	# print('X.shape: ', X.shape)
	# print('y.shape: ', y.shape)

	# Apply the weight for total energy minimum points
	weights = np.array([])
	
	for i in scan_fit.index.tolist():
		if i > 0 and i < len(scan)-1:
			if (scan.loc[i, 'E_qm - enqm_0 (kcal/mol)'] < scan.loc[i-1, 'E_qm - enqm_0 (kcal/mol)'] and 
				scan.loc[i, 'E_qm - enqm_0 (kcal/mol)'] < scan.loc[i+1, 'E_qm - enqm_0 (kcal/mol)']):
				weights = np.append(weights, weight_minimums)
			else:
				weights = np.append(weights, 1)
		else:
			weights = np.append(weights, 1)

	# Apply the linear regression model
	if method == 'least-square':		
		reg = LinearRegression(fit_intercept=False)
	elif method == 'ridge':
		reg = Ridge(alpha=lasso_alpha, max_iter=10000000, fit_intercept=False)
	elif method == 'lasso':
		reg = Lasso(alpha=lasso_alpha, max_iter=10000000, fit_intercept=False)
	elif method == 'lassocv':
		reg = LassoCV(eps=1e-5, max_iter=10000000, fit_intercept=False)
	else:
		print("The method is not specified or not implemented, using the default least square method.")
		reg = LinearRegression(fit_intercept=False)
	
	reg.fit(X, y, sample_weight=weights)

	if method == 'lassocv':
		print('Optimal penalty coefficient from Lasso CV regression: %s' % reg.alpha_)

	# Create a nested dictionary with the linear regression data
	lr_data = {}	# {i: {'tors': str, 'atoms': list, 'constants': list}}

	for i in range(int(len(reg.coef_)/6)):
		lr_data[i] = {}
		lr_data[i]["tors"] = X.columns.values.tolist()[6*i][:len(X.columns.values.tolist()[6*i])-3]
		lr_data[i]["atoms"] = []
		lr_data[i]["atoms"].extend([elem for elem in lr_data[i]["tors"].split('-')])
		lr_data[i]["constants"] = []
		for j in range(6):
			lr_data[i]["constants"].append(reg.coef_[6*i+j])

	# Write the linear regression coefficients in the topology file
	write_itp_file(topfile, lr_data)

	# Plot the results

	# Fitted dihedral energy
	y_plot = np.dot(X.values, reg.coef_)

	# Dihedral angles
	x_fit = scan_fit['Dihedral (rad)'].copy()
	x_full = scan['Dihedral (rad)'].copy()

	# X-axis points
	x_plot1 = np.linspace(x_fit.min(), x_fit.max(), 300)
	x_plot2 = np.linspace(x_full.min(), x_full.max(), 300)

	# Interpolate required angles
	smooth_dih = interp1d(x_fit, y_plot, kind='cubic')
	smooth_nben_fit = interp1d(x_fit, scan_fit['Non-bonded - nben0 (kcal/mol)'].copy(), kind='cubic')

	# Interpolate all angles available
	smooth_nben = interp1d(x_full, scan['Non-bonded - nben0 (kcal/mol)'].copy(), kind='cubic')
	smooth_qmE = interp1d(x_full, scan['E_qm - enqm_0 (kcal/mol)'].copy(), kind='cubic')
	smooth_qmDih = interp1d(x_full, scan['U_tors-qm (kcal/mol)'].copy(), kind='cubic')

	fig, ((ax1), (ax2)) = plt.subplots(ncols=1, nrows=2, figsize=(7, 5))

	# Total energy
	ax1.plot(x_plot2*180/np.pi, smooth_qmE(x_plot2), '--', c="tab:purple", label="Gaussian total energy")
	ax1.scatter(x_full*180/np.pi, scan['E_qm - enqm_0 (kcal/mol)'].copy(), c="tab:purple")

	ax1.plot(x_plot1*180/np.pi, smooth_dih(x_plot1) + smooth_nben_fit(x_plot1), '--', c="tab:green", label="Classical total energy")
	ax1.scatter(x_fit*180/np.pi, scan_fit['E_qm - enqm_0 (kcal/mol)'].copy(), c="tab:green", marker='x', s=50)

	# Dihedral energy
	ax2.plot(x_plot2*180/np.pi, smooth_qmDih(x_plot2), '--', c="tab:red", label="Gaussian torsional energy")
	ax2.scatter(x_full*180/np.pi, scan['U_tors-qm (kcal/mol)'].copy(), c="tab:red")

	ax2.plot(x_plot1*180/np.pi, smooth_dih(x_plot1), '--', c="tab:orange", label="Fitted torsional energy")
	ax2.scatter(x_fit*180/np.pi, y_plot, c="tab:orange", marker='x', s=50)

	# Plot setup
	ax1.set_ylabel('Energy (kcal/mol)')
	ax1.legend(frameon=False)
	# ax1.set_ylim(-0.3, 2.3)

	ax2.set_xlabel(r'Angle ($^o$)')
	ax2.set_ylabel('Energy (kcal/mol)')
	ax2.legend(frameon=False)

	plt.tight_layout()
	plt.savefig("linear_regression.png", bbox_inches='tight', format='png', dpi=600)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Performs a linear regression to determine torsional dihedral coefficients.")
	parser.add_argument("gaussianlogfile", help="the gaussian log file.")
	parser.add_argument("xyzrotationsfile", help="file with all torsional scan configurations.")
	parser.add_argument("topfile", help="the topology file.")
	parser.add_argument("txtfile", help="DICE .txt file")
	parser.add_argument("dfrfile", help="DICE .dfr file")
	parser.add_argument("a1", type=int, help="first atom defining the reference dihedral.")
	parser.add_argument("a2", type=int, help="second atom defining the reference dihedral.")
	parser.add_argument("a3", type=int, help="third atom defining the reference dihedral.")
	parser.add_argument("a4", type=int, help="fourth atom defining the reference dihedral.")
	parser.add_argument("npoints", type=int, help="number of configurations during the scan.")
	parser.add_argument("--method", "-m", help="the method employed in the linear regression (least-square, ridge, lasso, lassocv).", 
						default='least-square')
	parser.add_argument("--alpha", type=float, help="the coefficient multiplying L1 penalty in Lasso linear regression. Default is 0.01.", default=0.01)
	parser.add_argument("--weight", "-w", type=float, help="the weight given to total energy minima points. Default is 1.", default=1)
	parser.add_argument("--cutoff", "-c", type=float, help="minimum atomic distance tolerated. Default is 0.5.", default=0.5)
	parser.add_argument("--remove-overlap", "-r", help="remove the overlapping configurations.", action='store_true')
	parser.add_argument("--max-barrier", "-b", help="limit the torsional barriers to the provided value (in kcal/mol).", default=None)
	parser.add_argument("--fit-from-total", "-t", help="fit the torsional angle using the total energy.", action='store_true', default=False)
	parser.add_argument("--angles", nargs='+', type=float, help="all dihedrals (in degrees) to be used in the fit (include the max and min angles).")

	args = parser.parse_args()

	data = get_data(args.gaussianlogfile, args.txtfile, args.dfrfile, args.a1, args.a2, args.a3, args.a4)

	dih_fit, dih = get_dihedrals(args.xyzrotationsfile, args.topfile, args.a1, args.a2, args.a3, args.a4, args.npoints, args.angles)

	if args.remove_overlap:
		data, dih_fit = remove_overlap(data, dih_fit, args.topfile, args.xyzrotationsfile, args.cutoff, args.npoints)
	else:
		_, _ = check_overlap(args.topfile, args.xyzrotationsfile, args.cutoff, args.npoints)

	scan_fit, scan = rescale_data(data, args.angles, args.max_barrier)

	# linear_regression(scan, scan_fit, dih, dih_fit, args.topfile, args.method, args.alpha, args.weight, args.fit_from_total)
	linear_regression(scan, scan_fit, dih_fit, args.topfile, args.method, args.alpha, args.weight, args.fit_from_total)
