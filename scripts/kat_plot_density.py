#!/usr/bin/env python3

import argparse
import numpy as np
import scipy.ndimage as ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import kat_plot_colormaps as cmaps

from kat_plot_misc import *

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
                    help="The plot file type to create (default is based on " \
                    "given output name).")
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
parser.add_argument("--contours", choices=["none", "normal", "smooth"],
                    default="normal")
parser.add_argument("--not_rasterised", dest="rasterised",
                    action="store_false",
                    help="Don't rasterise graphics (slower).")
parser.add_argument("--dpi", type=int, default=300,
                    help="Resolution in dots per inch of output graphic.")
parser.set_defaults(rasterised=True)
parser.add_argument("-v", "--verbose", dest="verbose",
                    action="store_true",
                    help="Print extra information")
parser.set_defaults(verbose=False)

args = parser.parse_args()
# ----- end command line parsing -----

# load header information
input_file = open(args.matrix_file)

header = readheader(input_file)

if args.title is not None:
    title = args.title
elif "Title" in header:
    title = header["Title"]
else:
    title = "Density Plot"

if args.x_label is not None:
    x_label = args.x_label
elif "XLabel" in header:
    x_label = header["XLabel"]
else:
    x_label = "X"

if args.y_label is not None:
    y_label = args.y_label
elif "YLabel" in header:
    y_label = header["YLabel"]
else:
    y_label = "Y"

if args.z_label is not None:
    z_label = args.z_label
elif "ZLabel" in header:
    z_label = header["ZLabel"]
else:
    z_label = "Z"

matrix = np.loadtxt(input_file)
input_file.close()
if args.verbose:
    print("{:d} by {:d} matrix file loaded.".format(matrix.shape[0],
                                                    matrix.shape[1]))

if args.contours == "smooth":
    matrix_smooth = ndimage.gaussian_filter(matrix, sigma=2.0, order=0)

if args.x_max is None or args.y_max is None or args.z_max is None:
    # find peaks
    msum = np.sum(matrix)
    xsums = np.sum(matrix, 0)
    ysums = np.sum(matrix, 1)
    peakx = findpeaks(xsums)
    peaky = findpeaks(ysums)
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
    for i in range(1, len(xsums), int(len(xsums)/40) + 1):
        if np.sum(xsums[:i]) >= msum * 0.995:
            xmax = i
            break
    for i in range(1, len(ysums), int(len(ysums)/40) + 1):
        if np.sum(ysums[:i]) >= msum * 0.995:
            ymax = i
            break

    zmax = np.max(peakz) * 1.1

    if args.verbose:
        print("Automatically detected axis limits:")
        print("xmax: ", xmax)
        print("ymax: ", ymax)
        print("zmax: ", zmax)

if args.x_max is not None:
    xmax = args.x_max
if args.y_max is not None:
    ymax = args.y_max
if args.z_max is not None:
    zmax = args.z_max

plt.figure(num = None, figsize=(args.width, args.height))

pcol = plt.pcolormesh(matrix, vmin=0, vmax=zmax, cmap=cmaps.viridis,
                      rasterized=args.rasterised)
plt.axis([0,xmax,0,ymax])
cbar = plt.colorbar()
cbar.set_label(z_label)
cbar.solids.set_rasterized(args.rasterised)
levels = np.arange(zmax/8, zmax, zmax/8)
if args.contours == "normal":
    plt.contour(matrix, colors="white", alpha=0.6, levels=levels)
elif args.contours == "smooth":
    plt.contour(matrix_smooth, colors="white", alpha=0.6, levels=levels)

plt.title(title)
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.grid(True, color="white", alpha=0.2)

if args.output_type is not None:
    output_name = args.output + '.' + args.output_type
else:
    output_name = args.output

plt.savefig(correct_filename(output_name), dpi=args.dpi)

