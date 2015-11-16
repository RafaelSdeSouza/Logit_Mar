#!/usr/bin/env python

import pickle
import corner
import argparse
from os.path import splitext

# Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filename', help='data file name',
                    default='logit.p')
args = parser.parse_args()

# Read stan samples
with open(args.filename, 'rb') as f:
    samples = pickle.load(f)

ndims = samples.shape
truths = None

figure = corner.corner(samples[:,:,:].reshape(ndims[0] * ndims[1], ndims[2]),
                       quantiles=[0.16, 0.5, 0.84],
                       labels=[r"$\alpha$"] +
                       [r"$\beta_{0}$".format(i) for i in range(ndims[2]-1)],
                       truths=truths,
                       show_titles=True)
figure.savefig("{0}.png".format(splitext(args.filename)[0]))
