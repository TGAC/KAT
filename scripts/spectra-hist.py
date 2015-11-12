#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

from findpeaks import *
from header import *

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Creates K-mer spectra plot from one of more histograms.")

parser.add_argument("histo_files", type=str, nargs='+',
                    help="The input histogram file from KAT")

parser.add_argument("-o", "--output", type=str, default="kat-spectra-hist",
                    help="The path to the output file.")
parser.add_argument("-p", "--output_type", type=str,
                    help="The plot file type to create (default is based on given output name).")
parser.add_argument("-t", "--title", type=str,
                    help="Title for plot")
parser.add_argument("-a", "--x_label", type=str,
                    help="Label for x-axis")
parser.add_argument("-b", "--y_label", type=str,
                    help="Label for y-axis")
parser.add_argument("-L", "--legend_labels", type=str,
                    help="Comma separated list of labels for legend")
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
parser.add_argument("-m", "--x_logscale", dest="x_logscale", action="store_true",
                    help="X-axis is logscale. Overrides x_min and x_max")
parser.set_defaults(x_logscale=False)
parser.add_argument("-n", "--y_logscale", dest="y_logscale", action="store_true",
                    help="Y-axis is logscale. Overrides y_min and y_max")
parser.set_defaults(y_logscale=False)
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                    help="Print extra information")
parser.set_defaults(verbose=False)

args = parser.parse_args()
# ----- end command line parsing -----

if args.verbose:
    print "Input files given: {:d}".format(len(args.histo_files))

headers = []
x = []
y = []
for histo_file in args.histo_files:
    input_file = open(histo_file)

    header = readheader(input_file)

    matrix = np.loadtxt(input_file)
    input_file.close()
    if args.verbose:
        print "{:d} by {:d} matrix file loaded.".format(matrix.shape[0],
                                                        matrix.shape[1])
    headers.append(header)
    x.append(matrix[:,0])
    y.append(matrix[:,1])

if args.title is not None:
    title = args.title
elif "Title" in header:
    title = headers[0]["Title"]
else:
    title = "Spectra Copy Number Plot"

if args.x_label is not None:
    x_label = args.x_label
elif "XLabel" in header:
    x_label = headers[0]["XLabel"]
else:
    x_label = "X"

if args.y_label is not None:
    y_label = args.y_label
elif "YLabel" in header:
    y_label = headers[0]["YLabel"]
else:
    y_label = "Y"

# find limits
if args.x_max is None or args.y_max is None:
    xmax = map(len, x)
    ysum = map(np.sum, y)
    ymax = map(np.max, y)
    for i in range(len(x)):
        peakx = findpeaks(y[i])
        peakx = peakx[peakx != 1]
        peaky = y[i][peakx]

        for j in range(1, xmax[i], xmax[i]/1000 + 1):
            if np.sum(y[i][:j]) >= ysum[i] * 0.999:
                xmax[i] = j
                break

        ymax[i] = np.max(peaky) * 1.1

    xmax = max(xmax)
    ymax = max(ymax)
    if args.verbose:
        print "Automatically detected axis limits:"
        print "xmax: ", xmax
        print "ymax: ", ymax

if args.x_max is not None:
    xmax = args.x_max
if args.y_max is not None:
    ymax = args.y_max

plt.figure(num = None, figsize=(args.width, args.height))

legend_labels = []
if args.legend_labels is not None:
    legend_labels = args.legend_labels.split(",")
if len(legend_labels) >= len(x):
    labels = legend_labels
else:
    labels = map(lambda s: s.split("/")[-1], args.histo_files)

colours = ["#cc0000",
           "#75507b",
           "#3465a4",
           "#73d216",
           "#c17d11",
           "#f57900",
           "#edd400"]
for xt,yt,lb,i in zip(x,y,labels,range(len(x))):
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

plt.savefig(output_name, dpi=300)
