#!/usr/bin/python

import pandas as pd
import pystan
import argparse
from os.path import abspath, join, basename

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
args = parser.parse_args()

# read data into pandas data frame
datadir = abspath('../data')
datafile = join(datadir, "matched.txt")
df = pd.read_csv(datafile, sep=" *", index_col=False)

print(df.columns)

# Stan fit
with open(args.filename) as f:
    fit = pystan.stan(model_code=f.read(), data=data, iter=args.niter,
        chains=3, n_jobs=100)

a = fit.extract(permuted=False)

# pickle stan fit object
with open(join(
    datadir, "{0}.p".format(splitext(basename(args.filename))[0])), "wb") as f:
    pickle.dump(a, f, pickle.HIGHEST_PROTOCOL)
