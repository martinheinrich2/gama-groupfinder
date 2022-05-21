#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find Galaxy groups from galaxy links
"""

from timeit import default_timer as timer
import numpy as np
import pandas as pd
import os


def get_lastpair(gal_array):
    """ Function to get the last galaxy pair in the array

    Parameter:
        gal_array (numpy array): [array of all remaining galaxies]

    Returns:
        [numpy array]: [p_arr is last galaxy pair in the array]
    """
    cataid_1 = gal_array[-1][0]
    cataid_2 = gal_array[-1][1]
    p_index = np.where(gal_array == cataid_1)[0][0]
    p_arr = np.column_stack((p_index, cataid_1))
    p_arr = np.column_stack((p_arr, cataid_2))
    return p_arr


def get_neighbours(data_array, test_array):
    """"
    Function to return array of neighbours from work_data and array
    of galaxies

    Parameter: data_array and test_array

    Returns: array of indices from neighbour galaxies
    """
    new_index = np.empty([0, 1])
    new_array = np.empty([0, 2])
    n1_boolean = np.isin(data_array, test_array[0])
    n2_boolean = np.isin(data_array, test_array[1])
    n_1 = np.array(np.where(n1_boolean))
    n_2 = np.array(np.where(n2_boolean))
    if (len(n_1[0]) or len(n_2[0])) != 0:
        nnew_1 = n_1[0]
        nnew_2 = n_2[0]
        new_index = np.concatenate((nnew_1, nnew_2))
        new_array = data_array[new_index]
        new_index = np.unique(new_index)
        new_index - np.sort(new_index)
        new_array = np.unique(new_array)
        new_array = np.sort(new_array)
    return new_index, new_array


start1 = timer()
path = 'results/galaxy-links/'
# get list of all files in folder
filelist = [file for file in os.listdir(path) if file.endswith(".csv")]
for infile in filelist:
    start = timer()
    # split filename into components and create filename for results file
    f1 = infile.split("_")
    f2 = str("freq_") + f1[1] + str("_") + f1[2]
    export_filenameR = ('./results/galaxy-richness/' + f2)
    export_filenameG = ('./results/galaxy-groups/' + f2)
    print("Current file is: " + export_filenameR)
    my_data = np.genfromtxt(path + infile, delimiter=',', skip_header=1)
    my_data = np.delete(my_data, 0, 1)
    # sort data and remove duplicates
    foo1 = np.sort(my_data, axis=1)
    foo2 = np.unique(foo1, axis=0)
    # create empty arrays with two columns for CATAID and group
    groups_1 = np.empty([0, 2])
    groups_2 = np.empty([0, 2])
    work_data = np.copy(foo2.astype(int))
    new_arr = np.empty([0, 2])
    gnew_arr = np.empty([0, 2])
    gnew_ind = np.empty([0, 1])
    g_index = np.empty([0, 1])
    g_count = 1
    # start group finder loop
    while len(work_data) > 0:
        # get row number and CATAID of last row
        last_row = get_lastpair(work_data)
        g_array = last_row[0][1:, ]
        g_length = len(g_array)
        g_index = np.array(last_row[0][0])
        g_index = np.reshape(g_index, (1,))
        # delete last row from work array to avoid double detection
        work_data = np.delete(work_data, -1, 0)
        # get first array of index and galaxy IDs of all neighbours
        gnew_ind, gnew_arr = get_neighbours(work_data, g_array)
        # now find all neighbours
        while len(gnew_arr) > 0:
            g_array = np.append(g_array, gnew_arr, axis=0)
            g_index = np.append(g_index, gnew_ind, axis=0)
            # make arrays unique
            g_array = np.unique(g_array, axis=0)
            g_index = np.unique(g_index, axis=0)
            # delete rows of neighbours from array
            work_data = np.delete(work_data, gnew_ind, axis=0)
            gnew_ind, gnew_arr = get_neighbours(work_data, g_array)
            g_length = len(g_array)
        # assign group number and add galaxies to final groups array
        groups_1 = np.full((g_length,), g_count)
        groups_1 = np.column_stack((groups_1, g_array))
        groups_2 = np.append(groups_2, groups_1, axis=0)
        g_count += 1
        # print("\r", len(work_data), end="", flush=True)
    final_groups = groups_2.astype('int32')
    # create dataframe with all galaxy groups and add column names
    groups = pd.DataFrame(final_groups)
    groups.columns = ['group', 'CATAID']
    # create dataframe with group sizes
    g_freq = groups.groupby('group').sum()
    g_freq = groups['group'].value_counts(dropna=False)
    g_freq = pd.DataFrame(groups.group.value_counts(dropna=False))
    # create dataframe with richness of groups (group frequency)
    freq = pd.DataFrame(g_freq.group.value_counts(dropna=False))
    freq['richness'] = freq.index
    new_freq = freq.rename(columns={'group': 'count'})
    pd.DataFrame.to_csv(new_freq, export_filenameR, index=True)
    pd.DataFrame.to_csv(groups, export_filenameG, index=True)
    stop = timer()
    print("\n", stop - start)
    groups = groups.reset_index(drop=True)

stop1 = timer()
print("\n", stop1 - start1)
