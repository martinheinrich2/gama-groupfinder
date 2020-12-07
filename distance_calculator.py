# -*- coding: utf-8 -*-

"""distance_calculator

This script reads a sample.csv file from the GAMA-Survey and calculates
the co-moving distances, logarithmic shift (zeta) and converts
coordinate degrees to radians. It exports two versions with different
amount of columns as CSV-files. The group finding process had better
performance with the reduced data.
The volume-limited sample.csv for the project was created using TOPCAT.
"""

import numpy as np
import pandas as pd
from astropy.cosmology import wCDM

"""
Set the cosmology with parameters for dark matter (Om0), dark energy (Ode0),
dark energy equation of state at all redshifts (cosmological constant with
value w0 = -1.0) and the hubble parameter at z = 0.
"""
cosmo = wCDM(H0=70, Om0=0.3, Ode0=0.7, w0=-1.0)

# Read data from CSV file into a dataframe
print('Reading volume limited sample ...')
file = r'data/volume_limited_sample_z03.csv'
my_data = pd.read_csv(file)

"""
Calculate the line-of-sight (radial) and transverse distance in MPC
using the functions comoving_distance from astropy.cosmology subpackage.
A new column is created with the results for all rows in the table.
"""

print('calculating comoving distances...')
my_data['CoDist'] = cosmo.comoving_distance(my_data['Z_CMB'])
print('calculating the comoving transverse distance...')
my_data['CoDistTran'] = cosmo.comoving_transverse_distance(my_data['Z_CMB'])

"""
z_cmb is the redshift of the observed galaxies in the CMB frame define c
as object in km/s and add new column for line of velocity v = c*ln(1+z_cmb)
this is actually v = c * zeta for logarithmic shift
zeta = ln(1+z) see also Baldry 2018a (https://arxiv.org/abs/1812.05135)
calculate the line of velocity
"""
print('calculating the line of velocity...')
c = 299792.458
my_data['line_v'] = (c * np.log(1 + my_data['Z_CMB']))

# calculate zeta = ln(1+z_cmb)
print('calculating zeta...')
my_data['zeta'] = (np.log(1 + my_data['Z_CMB']))

# convert RA and DEC from degree to radians and add new columns
print('converting degrees to radians...')
my_data['RA_rad'] = np.radians(my_data['RA'])
my_data['DEC_rad'] = np.radians(my_data['DEC'])

# collecting relevant columns in new table
reduced_data = my_data[['CATAID', 'RA', 'DEC', 'Z_CMB', 'CoDist', 'CoDistTran',
                        'line_v', 'zeta', 'RA_rad', 'DEC_rad']].copy()
print('Writing converted data to file...')
# uncomment line if more data columns are required
pd.DataFrame.to_csv(my_data, 'data/converted_sample.csv', index=False)

# Create a reduced data set with the required columns for the other scripts
# with columns:
# 0 CATAID, 1 RA, 2, DEC, 3 Z_CMB, 4 CoDist, 5 CoDistTran, 6 line_v,
# 7 zeta, 8 RA in radians, 9 DEC in radians
pd.DataFrame.to_csv(reduced_data, 'data/reduced_sample.csv', index=False)
