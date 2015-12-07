#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from kat_plot_misc import *

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Creates K-mer spectra plot from selected rows and/or " \
    "columns in a \"comp\" matrix.")

parser.add_argument("matrix_file", type=str,
                    help="The input matrix file from KAT")

parser.add_argument("-o", "--output", type=str, default="kat-spectra-mx",
                    help="The path to the output file.")
parser.add_argument("-p", "--output_type", type=str,
                    help="The plot file type to create (default is based on " \
                    "given output name).")
parser.add_argument("-t", "--title", type=str, default="Spectra MX Plot",
                    help="Title for plot")
parser.add_argument("-a", "--x_label", type=str,
                    help="Label for x-axis")
parser.add_argument("-b", "--y_label", type=str,
                    help="Label for y-axis")
parser.add_argument("-r", "--x_min", type=int, default=0,
                    help="Minimum value for x-axis")
parser.add_argument("-s", "--y_min", type=int, default=0,
                    help="Minimum value for y-axis")
parser.add_argument("-x", "--x_max", type=int,
                    help="Maximum value for x-axis")
parser.add_argument("-y", "--y_max", type=int,
                    help="Maximum value for y-axis")
parser.add_argument("-w", "--width", type=int, default=8,
                    help="Width of canvas")
parser.add_argument("-l", "--height", type=int, default=6,
                    help="Height of canvas")
parser.add_argument("-i", "--intersection", dest="intersection",
                    action="store_true",
                    help="Activate intersection mode, which plots the " \
                    "shared and exclusive content found in the matrix.")
parser.set_defaults(intersection=False)
parser.add_argument("-c", "--list", type=str,
                    help="The list of columns or rows to select from the " \
                    "matrix (overrides -i)")
parser.add_argument("-e", "--exc_cutoff_d1", type=int, default=1,
                    help="If in intersection mode, the level at which " \
                    "content for dataset 1 is considered exclusive or shared")
parser.add_argument("-f", "--exc_cutoff_d2", type=int, default=1,
                    help="If in intersection mode, the level at which " \
                    "content for dataset 2 is considered exclusive or shared")
parser.add_argument("-m", "--x_logscale", dest="x_logscale",
                    action="store_true",
                    help="X-axis is logscale. Overrides x_min and x_max")
parser.set_defaults(x_logscale=False)
parser.add_argument("-n", "--y_logscale", dest="y_logscale",
                    action="store_true",
                    help="Y-axis is logscale. Overrides y_min and y_max")
parser.set_defaults(y_logscale=False)
parser.add_argument("--dpi", type=int, default=300,
                    help="Resolution in dots per inch of output graphic.")
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

matrix = np.loadtxt(input_file)
input_file.close()
if args.verbose:
    print("{:d} by {:d} matrix file loaded.".format(matrix.shape[0],
                                                    matrix.shape[1]))

# make rows/columns of data to plot
x = []
y = []
labels = []
if args.list is not None:
    if args.verbose:
        print("List given: {:s}".format(args.list))
    rowscols = []
    try:
        for str in args.list.split(','):
            if str[0] in "rc":
                rowscols.append((str[0], int(str[1:])))
            else:
                raise ValueError()
    except ValueError as e:
        sys.exit("Malformed string given as --list: " + args.list)
    for rowcol in rowscols:
        if rowcol[0] == 'r':
            y.append(matrix[rowcol[1],:])
            x.append(np.arange(len(matrix[rowcol[1],:])))
            labels.append("Row {:d}".format(rowcol[1]))
        elif rowcol[0] == 'c':
            y.append(matrix[:,rowcol[1]])
            x.append(np.arange(len(matrix[:,rowcol[1]])))
            labels.append("Column {:d}".format(rowcol[1]))
elif args.intersection:
    if args.verbose:
        print("Intersection mode.")
        print("Dataset 1 cutoff: {:d}".format(args.exc_cutoff_d1))
        print("Dataset 2 cutoff: {:d}".format(args.exc_cutoff_d2))
    y_exc_d1 = np.transpose(np.sum(matrix[:,:args.exc_cutoff_d1],1))
    x_exc_d1 = np.arange(len(y_exc_d1))
    y_sha_d1 = np.transpose(np.sum(matrix[args.exc_cutoff_d2:,
                                          args.exc_cutoff_d1:],1))
    x_sha_d1 = np.arange(args.exc_cutoff_d2, len(y_exc_d1))
    y_exc_d2 = np.sum(matrix[:args.exc_cutoff_d2,:],0)
    x_exc_d2 = np.arange(len(y_exc_d2))
    y_sha_d2 = np.sum(matrix[args.exc_cutoff_d2:,
                             args.exc_cutoff_d1:],0)
    x_sha_d2 = np.arange(args.exc_cutoff_d1, len(y_exc_d2))
    x = [x_exc_d1, x_sha_d1, x_exc_d2, x_sha_d2]
    y = [y_exc_d1, y_sha_d1, y_exc_d2, y_sha_d2]
    labels = ["Dataset 1 exclusive content", "Dataset 1 shared content",
              "Dataset 2 exclusive content", "Dataset 2 shared content"]
else:
    sys.exit("Error: Either --list of --intersection must be given.")

# find limits
if args.x_max is None or args.y_max is None:
    xmax = list(map(len, x))
    ysum = list(map(np.sum, y))
    ymax = list(map(np.max, y))
    for i in range(len(x)):
        peakx = findpeaks(y[i])
        peakx = peakx[peakx != 1]
        peaky = y[i][peakx]

        for j in range(1, xmax[i], int(xmax[i]/1000) + 1):
            if np.sum(y[i][:j]) >= ysum[i] * 0.999:
                xmax[i] = j
                break

        ymax[i] = np.max(peaky) * 1.1

    xmax = max(xmax)
    ymax = max(ymax)
    if args.verbose:
        print("Automatically detected axis limits:")
        print("xmax: ", xmax)
        print("ymax: ", ymax)

if args.x_max is not None:
    xmax = args.x_max
if args.y_max is not None:
    ymax = args.y_max

plt.figure(num = None, figsize=(args.width, args.height))

colours = ["#cc0000",
           "#75507b",
           "#3465a4",
           "#73d216",
           "#c17d11",
           "#f57900",
           "#edd400"]
for xt,yt,lb,i in zip(x,y,labels,list(range(len(x)))):
    plt.plot(xt, yt, label=lb, color=colours[i%len(colours)])

if args.x_logscale:
    plt.xscale("log")
if args.y_logscale:
    plt.yscale("log")

plt.axis([args.x_min, xmax, args.y_min, ymax])

plt.title(title)
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.grid(True, color="black", alpha=0.2)
if len(x) > 1:
    plt.legend(loc=1)

if args.output_type is not None:
    output_name = args.output + '.' + args.output_type
else:
    output_name = args.output

plt.savefig(correct_filename(output_name), dpi=args.dpi)
