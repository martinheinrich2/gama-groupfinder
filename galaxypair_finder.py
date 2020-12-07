#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
find galaxy pairs within given line-of-sight and transverse linking lengths
and write link-number, CATAID1, CATAID2 in file
data structure in numpy array is:
0 CATAID, 1 RA, 2, DEC, 3 Z_CMB, 4 CoDist, 5 CoDistTran, 6 line_v,
7 zeta, 8 RA in radians, 9 DEC in radians

Uses numba module to make the code faster
"""
import numpy as np
import pandas as pd
from timeit import default_timer as timer
from numba import jit


@jit(nopython=True, parallel=True)
def calc_seplos(line_v1, line_v2):
    """Function to calculate the line-of-sight separation between galaxies

    Parameters: line-of-sight velocity of galaxy pair, line_v1 and
    line_v2

    Returns: sep_los, the separation of a galaxy pair in km/s
    """
    sep_los = np.abs(line_v1 - line_v2)
    return sep_los


@jit(nopython=True, parallel=True)
def calc_skysep(dec1, dec2, ra1, ra2):
    """Function to calculate the sky separation between two galaxies as
    sky_sep = np.arccos(((np.sin(DEC1)*np.sin(DEC2))+np.cos(DEC1)
    *np.cos(DEC2))*np.cos(RA1-RA2))
    from their DEC and RA coordinates in radians

    Parameters: dec1, dec2, ra1, ra2

    Returns: sky_sep, sky separation in radians
    """
    sky_sep = np.arccos(np.minimum(1, ((np.sin(dec1) * np.sin(dec2))
                                       + np.cos(dec1) * np.cos(dec2))
                                   * np.cos(ra1 - ra2)))
    return sky_sep


@jit(nopython=True, parallel=True)
def con_check(data, x_1, los, trans):
    """
    Function to add separations and linking conditions, looks if
    there a galaxy has neighbours.

    Parameters: data, dataframe of galaxies
                x_1, rownumber in dataframe
                los, line-of-sight velocity in km/s
                trans, transverse linking length in Mpc

    Returns:     A dataframe of galaxy CATAID's within linking lengths
    """
    # Add column of line-of-sight separation from function CalcSeplos
    data_seplos = (calc_seplos(data[x_1, 6], data[0:, 6]))
    """
    Add column of transverse separation from function CalcSkySep
    separation = theta * (CoDistTran[x] + CoDistTran)/2
    theta <- acos(sin(DEC1)sin(DEC2)+cos(DEC1)cos(DEC2)cos(RA1-RA2))
    first DEC then RA! see also wikipedia
    theta <- CalcSkySep(tmpData$RA_rad[x],tmpData$RA_rad,tmpData$DEC_rad[x],
    tmpData$DEC_rad) (wrong order of RA and DEC)
    """
    data_septrans = (calc_skysep(data[x_1, 9], data[0:, 9], data[x_1, 8],
                                 data[0:, 8]) * (data[x_1, 5] + data[0:, 5]) / 2)
    newdata = np.column_stack((data, data_seplos))
    newdata = np.column_stack((newdata, data_septrans))
    # Add column to check if both conditions are satisfied
    # (if only one true result, then there are no neighbours)
    tfdata = (newdata[0:, 10] <= los) & (newdata[0:, 11] <= trans)
    newdata = np.column_stack((newdata, tfdata))
    # filter only true conditions
    tmp_data2 = np.where(newdata[0:, 12])
    # filter only column CATAID
    tmp_data3 = newdata[0:, 0][tmp_data2]
    # return numpy array with CATAID values
    return tmp_data3


@jit(nopython=True, parallel=True)
def get_galaxypairs(los_1, trans_1, row_count, result_arr, tmp_data):
    """loop over all rows in dataframe and create and return dataframe 
    with galaxy pairs, CATAID1, CATAID2 and link-number"""
    for i in range(row_count):
        cataid_1 = tmp_data[i, 0]
        tmp_data1 = con_check(tmp_data, i, los_1, trans_1)
        if tmp_data1.size > 1:
            # delete element where CATAID1 = CATAID, e.g. a pair
            # of (12345 = 12345)
            tmp_data1 = tmp_data1[tmp_data1 != cataid_1]
            # get length of array and create arrays of equal lengths
            array_length = len(tmp_data1)
            CATAID1_arr = np.full((array_length,), cataid_1)
            # combine all arrays to one
            new_arr = np.column_stack((CATAID1_arr, tmp_data1))
            # create final array
            result_arr = np.append(result_arr, new_arr, axis=0)
        else:
            continue
    return result_arr


# define variable for the linking length in line-of-sigth (in km/s)
ll_los = float(500)
# load data into numpy array
my_data = np.genfromtxt('data/reduced_sample.csv', delimiter=',', skip_header=1)
# adjust range of transverse linking length, e.g. 0.1 to 2.0 Mpc in
# 0.1 Mpc steps
for l_trans in range(1, 21, 1):
    ll_trans = l_trans / 10
    print("line-of-sight", ll_los, "transverse =", ll_trans)
    # from timeit import default_timer as timer
    start = timer()
    tmp_data = my_data.copy()
    # create empty array of two coulumns for results
    final_arr = np.empty([0, 2])
    # create dataframe for resulting group pairs
    link_df = pd.DataFrame(columns=['CATAID1', 'CATAID2'])
    # get the number of rows in the dataframe
    rowcount = tmp_data.shape[0]
    # get galaxy pairs from function
    final_pairs = get_galaxypairs(ll_los, ll_trans, rowcount, final_arr, tmp_data)
    link_df = pd.DataFrame(data=final_pairs)
    link_df = link_df.astype(int)
    link_df.columns = ['CATAID1', 'CATAID2']
    # create filename from linking lengths
    filename = str('links_' + str(int(ll_los)) + '_' + str(ll_trans) + '.csv')
    export_filename = ('./results/galaxy-links/' + filename)
    # write data to file
    pd.DataFrame.to_csv(link_df, export_filename, index=True)
    end = timer()
    print("\n Execution time: ", end - start)
print("\n processed all linking lengths")
