#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from kat_plot_misc import *

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
parser.add_argument("-c", "--coverage_list", type=str,
                    help="Comma separated string listing coverage levels " \
                    "to show in plot (overrides -i)")
parser.add_argument("-u", "--cumulative", dest="cumulative",
                    action="store_true",
                    help="Plot cumulative distribution of kmers")
parser.set_defaults(cumulative=False)
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
    title = "k-mer comparison plot"

if args.x_label is not None:
    x_label = args.x_label
else:
    x_label = "k-mer multiplicity"

if args.y_label is not None:
    y_label = args.y_label
else:
    y_label = "Number of distinct k-mers"

matrix = np.loadtxt(input_file)
if header["Transpose"] == '1':
    matrix = np.transpose(matrix)
input_file.close()
if args.verbose:
    print("{:d} by {:d} matrix file loaded.".format(matrix.shape[0],
                                                    matrix.shape[1]))

mincov = 1 if args.ignore_absent else 0
covbands = args.max_dup

colours = ["#000000",
           "#ef2929",
           "#ad7fa8",
           "#8ae234",
           "#729fcf",
           "#f2c27e",
           "#fcaf3e",
           "#fce94f"]
if mincov > 0:
    colours = colours[1:]
    xamount = 0.65
else:
    xamount = 0.99

# leave only coverage levels we are interested in
last_row = np.matrix(np.sum(matrix[(mincov+covbands):,:], 0))
matrix = np.concatenate([matrix[mincov:(mincov+covbands),:], last_row], 0)

# find limits
if args.x_max is None or args.y_max is None:
    totals = np.squeeze(np.asarray(np.sum(matrix, 0)))
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

matrix = matrix[:,:xmax]

plt.figure(num = None, figsize=(args.width, args.height))
plt.axis([0,xmax,0,ymax])
x = list(range(xmax))
labels = ["{:d}x".format(l) for l in range(mincov, mincov+covbands+1)]
labels[-1] = "{:s}+".format(labels[-1])
bar = plt.bar(x, np.squeeze(np.asarray(matrix[0,:])),
              color=colours[0],
              linewidth=0.1,
              edgecolor=colours[0],
              width=1,
              label=labels[0])
for level in range(1, covbands+1):
    bar = plt.bar(x, np.squeeze(np.asarray(matrix[level,:])),
                  bottom=np.squeeze(np.asarray(np.sum(matrix[:level,:], 0))),
                  color=colours[level%len(colours)],
                  linewidth=0.1,
                  edgecolor=colours[level%len(colours)],
                  width=1,
                  label=labels[level])

plt.title(wrap(title))
plt.xlabel(wrap(x_label))
plt.ylabel(wrap(y_label))
plt.grid(True, color="black", alpha=0.2)
plt.legend(loc=1)
plt.tight_layout()

if args.output_type is not None:
    output_name = args.output + '.' + args.output_type
else:
    output_name = args.output

plt.savefig(correct_filename(output_name), dpi=args.dpi)
