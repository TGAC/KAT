#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Creates a stacked histogram showing the level of duplication in an assembly.")

parser.add_argument("matrix_file", type=str,
                    help="The input matrix file from KAT")

parser.add_argument("-p", "--output_type", type=str, default="png",
                    choices=["png", "pdf", "eps"],
                    help="The plot file type to create.")
parser.add_argument("-o", "--output", type=str, default="kat-spectra-cn",
                    help="The path to the output file.")
parser.add_argument("-t", "--title", type=str, default="Spectra Copy Number Plot",
                    help="Title for plot")
parser.add_argument("-a", "--x_label", type=str, default="X",
                    help="Label for x-axis")
parser.add_argument("-b", "--y_label", type=str, default="Y",
                    help="Label for y-axis")
parser.add_argument("-x", "--x_max", type=int, default=1000,
                    help="Maximum value for x-axis")
parser.add_argument("-y", "--y_max", type=int, default=1000,
                    help="Maximum value for y-axis")
parser.add_argument("-w", "--width", type=int, default=8,
                    help="Width of canvas")
parser.add_argument("-l", "--height", type=int, default=6,
                    help="Height of canvas")
parser.add_argument("-i", "--ignore_absent", dest="ignore_absent", action="store_true",
                    help="Ignore K-mers in reads but absent from the assembly")
parser.set_defaults(ignore_absent=False)
parser.add_argument("-m", "--max_dup", type=int, default=6,
                    help="Maximum duplication level to show in plots")
parser.add_argument("-c", "--columns", type=str,
                    help="Comma separated string listing columns to show in plot (overrides -a)")
parser.add_argument("-u", "--cumulative", dest="cumulative", action="store_true",
                    help="Plot cumulative distribution of kmers")
parser.set_defaults(cumulative=False)
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
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
    print "{:d} by {:d} matrix file loaded.".format(matrix.shape[0],
                                                    matrix.shape[1])

