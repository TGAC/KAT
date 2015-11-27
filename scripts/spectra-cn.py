#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

from findpeaks import *
from header import *

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Creates a stacked histogram showing the level of " \
    "duplication in an assembly.")

parser.add_argument("matrix_file", type=str,
                    help="The input matrix file from KAT")

parser.add_argument("-o", "--output", type=str, default="kat-spectra-cn",
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
parser.add_argument("-x", "--x_max", type=int,
                    help="Maximum value for x-axis")
parser.add_argument("-y", "--y_max", type=int,
                    help="Maximum value for y-axis")
parser.add_argument("-w", "--width", type=int, default=8,
                    help="Width of canvas")
parser.add_argument("-l", "--height", type=int, default=6,
                    help="Height of canvas")
parser.add_argument("-i", "--ignore_absent", dest="ignore_absent",
                    action="store_true",
                    help="Ignore K-mers in reads but absent from the " \
                    "assembly")
parser.set_defaults(ignore_absent=False)
parser.add_argument("-m", "--max_dup", type=int, default=6,
                    help="Maximum duplication level to show in plots")
parser.add_argument("-c", "--columns", type=str,
                    help="Comma separated string listing columns to " \
                    "show in plot (overrides -a)")
parser.add_argument("-u", "--cumulative", dest="cumulative",
                    action="store_true",
                    help="Plot cumulative distribution of kmers")
parser.set_defaults(cumulative=False)
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

mincov = 1 if args.ignore_absent else 0
covbands = args.max_dup

colours = ["#888a85",
           "#ef2929",
           "#ad7fa8",
           "#729fcf",
           "#8ae234",
           "#e9b96e",
           "#fcaf3e",
           "#fce94f"]
if mincov > 0:
    colours = colours[1:]
    xamount = 0.65
else:
    xamount = 0.99

# leave only coverage levels we are interested in
last_column = np.transpose(np.matrix(np.sum(matrix[:,(mincov+covbands):], 1)))
matrix = np.concatenate([matrix[:,mincov:(mincov+covbands)], last_column], 1)

# find limits
if args.x_max is None or args.y_max is None:
    totals = np.squeeze(np.asarray(np.sum(matrix, 1)))
    xmax = len(totals) - 1
    ysum = np.sum(totals)
    ymax = np.max(totals)

    peakx = findpeaks(totals)
    peakx = peakx[peakx != 1]
    peaky = totals[peakx]

    for i in range(1, xmax, int(xmax/100) + 1):
        if np.sum(totals[:i]) >= ysum * xamount:
            xmax = i
            break
    ymax = np.max(peaky) * 1.1
    if args.verbose:
        print("Automatically detected axis limits:")
        print("xmax: ", xmax)
        print("ymax: ", ymax)

if args.x_max is not None:
    xmax = args.x_max
if args.y_max is not None:
    ymax = args.y_max

matrix = matrix[:xmax,:]

plt.figure(num = None, figsize=(args.width, args.height))
plt.axis([0,xmax,0,ymax])
x = list(range(xmax))
labels = ["{:d}x".format(l) for l in range(mincov, mincov+covbands+1)]
labels[-1] = "{:s}+".format(labels[-1])
bar = plt.bar(x, matrix[:,0],
              color=colours[0],
              linewidth=0.1,
              edgecolor=colours[0],
              width=1,
              label=labels[0])
for level in range(1, covbands+1):
    bar = plt.bar(x, matrix[:,level],
                  bottom=np.sum(matrix[:,:level], 1),
                  color=colours[level%len(colours)],
                  linewidth=0.1,
                  edgecolor=colours[level%len(colours)],
                  width=1,
                  label=labels[level])

plt.title(title)
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.grid(True, color="black", alpha=0.2)
plt.legend(loc=1)

if args.output_type is not None:
    output_name = args.output + '.' + args.output_type
else:
    output_name = args.output

plt.savefig(output_name, dpi=300)
