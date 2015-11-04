#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

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
parser.add_argument("-x", "--x_max", type=int, default=1000,
                    help="Maximum value for x-axis")
parser.add_argument("-y", "--y_max", type=int, default=1000,
                    help="Maximum value for y-axis")
parser.add_argument("-z", "--z_max", type=int, default=1000,
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

