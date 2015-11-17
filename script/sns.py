#!/usr/bin/env python

import pickle
import argparse
from os.path import splitext
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', help='data file name',
                    default='logit.p')
args = parser.parse_args()

# # Read stan samples
with open(args.filename, 'rb') as f:
    xvar, yvar, samples = pickle.load(f)

ndims = samples.shape
nx = len(xvar)

df = pd.DataFrame(samples[:,:,:nx+1].reshape(ndims[0] * ndims[1], nx+1),
        columns=[r"$\alpha$"] + [r"$\beta_{0}$".format(i) for i in range(nx)])

# g = sns.jointplot(x=r"$\alpha$", y=r"$\beta_0$", data=df, kind="kde")
g = sns.pairplot(df)
g.map_diag(sns.kdeplot)
g.map_offdiag(sns.kdeplot)
plt.savefig("{0}_sns.png".format(splitext(args.filename)[0]))

# print(df)
# sns.violinplot(y="color_gr", x="bpt", hue="bpt", data=df, split=True)
# plt.savefig("output.png")
