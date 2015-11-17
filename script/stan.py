#!/usr/bin/python

import pandas as pd
import numpy as np
import pystan
import argparse
from os.path import abspath, join, basename, splitext
import pickle

def normalize(df, var):
    """Normalize variables and uncertainties in place.

    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame.
    var : string
        Name of variable to be normalized.

    Returns
    -------
    None
    """
    std = np.std(df[var])
    df[var] = (df[var] - np.mean(df[var]))/std

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', help='data file name',
                    default='stan_models/logit.stan')
parser.add_argument('-n', '--niter', help='number of iterations',
                    type=int, default=1000)
parser.add_argument('-x', '--xvar', help='predictor variable(s)',
                    nargs='+', default=["lgm_tot_p50", "logM200_L"],
                    choices = ["lgm_tot_p50", "logM200_L", "RprojLW_Rvir",
                        "sfr_tot_p50", "color_gr"])
parser.add_argument('-c', '--chains', help='number of chains',
                    type=int, default=3)
args = parser.parse_args()

# read data into pandas data frame
datadir = abspath('../data')
datafile = join(datadir, "matched.txt")
df = pd.read_csv(datafile, sep=" *", index_col=False)


# Normalize variables
for v in args.xvar:
    normalize(df, v)

# build data dictionary for stan model
data = {}
data['n'] = len(df[args.xvar[0]])
data['k'] = len(args.xvar)
data['x'] = df[args.xvar].values
data['y'] = df["bpt"].values

# Stan fit
with open(args.filename) as f:
    fit = pystan.stan(model_code=f.read(), data=data, iter=args.niter,
                    chains=args.chains)

a = fit.extract(permuted=False)

# pickle stan fit object
with open(join(
    datadir, "{0}.p".format(splitext(basename(args.filename))[0])), "wb") as f:
    pickle.dump([args.xvar, "bpt", a], f, pickle.HIGHEST_PROTOCOL)
