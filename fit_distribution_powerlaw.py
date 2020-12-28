'''
Calculate powerlaw fit and parameters
Power laws are probability distributions with the form
p(x) proportional to x^-alpha
here we want to find some parameters and compare them later
'''
# %%
# powerlaw is not in anaconda, so do: conda install -c conda-forge powerlaw
import powerlaw
import numpy as np
import pandas as pd
import os

# set data path
path = 'results/galaxy-richness/'
# create empty dataframe
param_df = pd.DataFrame()

# get a list of all .csv files in the folder
filelist = [file for file in os.listdir(path) if file.endswith(".csv")]
# %%


def powerlawFitting(data):
    """
    Fits the data to a powerlaw, extracts and returns parameters

    Args:
            data (numpy array): numpy array of group richness, being the
            distribution of galaxy groups

    Returns:
            Powerlaw parameters, alpha [slope], sigma [standard error],
            x_min [data value beyond which the distribution was fitted],
            normalization [constant C], tailRatio [tail length/all elements]
    """
    # standard fit method is 'Likelihood' the maximum Likelihood estimation
    # it fits the parameters of the distribution to the data
    fit = powerlaw.Fit(data[:, 0])
    # get the number of values in array
    n = len(data[:, 0])
    #
    alpha = fit.alpha
    # standard error of the power law fit
    sigma = fit.sigma
    # the power law fit starts from the optimal xmin, if none provided, it will
    # be calculated
    x_min = fit.xmin
    # calculate normalisation C=(alpha-1)*xmin^(alpha-1) from the continuous case
    normalisation = (alpha - 1) * (x_min ** (alpha - 1))
    # minimal Kolmogorov-Smirnov distance between the data and the fit
    D = fit.D
    # length of tail and tail ratio
    n_tail = fit.n_tail
    tailRatio = n_tail / n

    # print(fit.distribution_compare('power_law', 'exponential'))
    return alpha, sigma, x_min, normalisation, tailRatio, n, D


# %%
for infile in filelist:
    f1 = infile.split("_")
    f2 = str(f1[2])
    f3 = f2.split(".")
    llos = f1[1]
    ltrans = f3[0] + str(".") + f3[1]
    print("Current file is: " + llos + " " + ltrans)
    data_array = np.genfromtxt(path + infile, delimiter=',', skip_header=1,
                               usecols=(1, 2))

    alpha, sigma, x_min, normalisation, tailRatio, n, D\
        = powerlawFitting(data_array)

    df_temp = pd.DataFrame([llos, ltrans, alpha, sigma, x_min,
                            normalisation, tailRatio, n, D])
    df_temp = df_temp.transpose()
    param_df = param_df.append(df_temp, ignore_index=True)
# %%
param_df.columns = ['llos', 'ltrans', 'alpha', 'sigma', 'xmin', 'normalisation',
                    'tailRatio', 'n', 'D']
param_df.to_csv('results/powerlaw_parameters.csv')
# %%
