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

parser.add_argument("-o", "--output", type=str, default="kat-density",
                    help="The path to the output file.")
parser.add_argument("-p", "--output_type", type=str,
                    help="The plot file type to create (default is based on given output name).")
parser.add_argument("-t", "--title", type=str,
                    help="Title for plot")
parser.add_argument("-a", "--x_label", type=str,
                    help="Label for x-axis")
parser.add_argument("-b", "--y_label", type=str,
                    help="Label for y-axis")
parser.add_argument("-c", "--z_label", type=str,
                    help="Label for z-axis")
parser.add_argument("-x", "--x_max", type=int,
                    help="Maximum value for x-axis")
parser.add_argument("-y", "--y_max", type=int,
                    help="Maximum value for y-axis")
parser.add_argument("-z", "--z_max", type=int,
                    help="Maximum value for z-axis")
parser.add_argument("-w", "--width", type=int, default=8,
                    help="Width of canvas")
parser.add_argument("-l", "--height", type=int, default=6,
                    help="Height of canvas")
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                    help="Print extra information")
parser.set_defaults(verbose=False)

args = parser.parse_args()
# ----- end command line parsing -----

# load header information
header_title = header_x = header_y = header_z = None
input_file = open(args.matrix_file)
for line in input_file:
    if line[0] == '#':
        if line[2:8] == "Title:":
            header_title = line[8:-1]
            if args.verbose:
                print header_title
        elif line[2:9] == "XLabel:":
            header_x = line[9:-1]
            if args.verbose:
                print header_x
        elif line[2:9] == "YLabel:":
            header_y = line[9:-1]
            if args.verbose:
                print header_y
        elif line[2:9] == "ZLabel:":
            header_z = line[9:-1]
            if args.verbose:
                print header_z
        elif line[0:-1] == "###":
            break
    else:
        break

if args.title is not None:
    title = args.title
elif header_title is not None:
    title = header_title
else:
    title = "Density Plot"

if args.x_label is not None:
    x_label = args.x_label
elif header_x is not None:
    x_label = header_x
else:
    x_label = "X"

if args.y_label is not None:
    y_label = args.y_label
elif header_y is not None:
    y_label = header_y
else:
    y_label = "Y"

if args.z_label is not None:
    z_label = args.z_label
elif header_z is not None:
    z_label = header_z
else:
    z_label = "Z"

matrix = np.loadtxt(input_file)
input_file.close()
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
    peaky = peaky[peaky != 1]
    peakz = matrix[peaky,:][:,peakx]

    # peakxv = xsums[peakx]
    # print "peakxv: ", peakxv
    # xmax = np.max(peakx[peakxv > (msum * 0.0005)]) * 2
    # peakyv = ysums[peaky]
    # print "peakyv: ", peakyv
    # ymax = np.max(peaky[peakyv > (msum * 0.0005)]) * 2

    xmax = len(xsums)
    ymax = len(ysums)
    for i in range(1, len(xsums), len(xsums)/40 + 1):
        if np.sum(xsums[:i]) >= msum * 0.995:
            xmax = i
            break
    for i in range(1, len(ysums), len(ysums)/40 + 1):
        if np.sum(ysums[:i]) >= msum * 0.995:
            ymax = i
            break

    zmax = np.max(peakz) * 1.1

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

plt.figure(num = None, figsize=(args.width, args.height))

pcol = plt.pcolormesh(matrix, vmin=0, vmax=zmax, cmap=cmaps.viridis, rasterized=True)
plt.axis([0,xmax,0,ymax])
cbar = plt.colorbar()
cbar.set_label(z_label)
levels = np.arange(zmax/4, zmax, zmax/8)
plt.contour(matrix_smooth, colors="white", alpha=0.6, levels=levels)

plt.title(title)
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.grid(True, color="white", alpha=0.2)

if args.output_type is not None:
    output_name = args.output + '.' + args.output_type
else:
    output_name = args.output

plt.savefig(output_name, dpi=300)
