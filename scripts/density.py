#!/usr/bin/env python

import argparse
import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import colormaps as cmaps

from findpeaks import *

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="""Create K-mer Density Plots.

Creates a scatter plot, where the density or "heat" at each point represents
the number of distinct K-mers at that point.  Typically this is used to
visualise a matrix produced by the "kat comp" tool to compare multiplicities
from two K-mer hashes produced by different NGS reads, or to visualise the GC
vs K-mer multiplicity matrices produced by the "kat gcp" tool.""")

parser.add_argument("matrix_file", type=str,
                    help="The input matrix file from KAT")

parser.add_argument("-p", "--output_type", type=str, default="png",
                    choices=["png", "pdf", "eps"],
                    help="The plot file type to create.")
parser.add_argument("-o", "--output", type=str, default="kat-density",
                    help="The path to the output file.")
parser.add_argument("-t", "--title", type=str, default="Density Plot",
                    help="Title for plot")
parser.add_argument("-a", "--x_label", type=str, default="X",
                    help="Label for x-axis")
parser.add_argument("-b", "--y_label", type=str, default="Y",
                    help="Label for y-axis")
parser.add_argument("-c", "--z_label", type=str, default="Z",
                    help="Label for z-axis")
parser.add_argument("-x", "--x_max", type=int,
                    help="Maximum value for x-axis")
parser.add_argument("-y", "--y_max", type=int,
                    help="Maximum value for y-axis")
parser.add_argument("-z", "--z_max", type=int,
                    help="Maximum value for z-axis")
parser.add_argument("-w", "--width", type=int, default=1024,
                    help="Width of canvas")
parser.add_argument("-l", "--height", type=int, default=1024,
                    help="Height of canvas")
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                    help="Print extra information")
parser.set_defaults(verbose=False)

args = parser.parse_args()
# ----- end command line parsing -----

matrix = np.loadtxt(args.matrix_file)
if args.verbose:
    print "{:d} by {:d} matrix file loaded.".format(matrix.shape[0],
                                                    matrix.shape[1])
matrix_smooth = ndimage.gaussian_filter(matrix, sigma=2.0, order=0)

if args.x_max is None or args.y_max is None or args.z_max is None:
    # find peaks
    msum = np.sum(matrix)
    xsums = np.sum(matrix, 0)
    ysums = np.sum(matrix, 1)
    peakx = findpeaks(xsums) + 1
    peaky = findpeaks(ysums) + 1
    # ignore peaks at 1
    peakx = peakx[peakx != 1]
    print peakx
    peaky = peaky[peaky != 1]
    print peaky
    peakz = matrix[peaky,:][:,peakx]

    # peakxv = xsums[peakx]
    # print "peakxv: ", peakxv
    # xmax = np.max(peakx[peakxv > (msum * 0.0005)]) * 2
    # peakyv = ysums[peaky]
    # print "peakyv: ", peakyv
    # ymax = np.max(peaky[peakyv > (msum * 0.0005)]) * 2

    xmax = len(xsums)
    ymax = len(ysums)
    for i in range(1, len(xsums), len(xsums)/20):
        if np.sum(xsums[:i]) >= msum * 0.995:
            xmax = i
            break
    for i in range(1, len(ysums), len(ysums)/20):
        if np.sum(ysums[:i]) >= msum * 0.995:
            ymax = i
            break


    zmax = np.max(peakz)

    if args.verbose:
        print "Automatically detected axis limits:"
        print "xmax: ", xmax
        print "ymax: ", ymax
        print "zmax: ", zmax

if args.x_max is not None:
    xmax = args.x_max
if args.y_max is not None:
    ymax = args.y_max
if args.z_max is not None:
    zmax = args.z_max

plt.pcolormesh(matrix, vmin=0, vmax=zmax, cmap=cmaps.viridis)
plt.axis([0,xmax,0,ymax])
plt.colorbar()
levels = np.arange(zmax/4, zmax, zmax/8)
plt.contour(matrix_smooth, colors="white", alpha=0.6, levels=levels)
plt.savefig(args.output + ".png")
