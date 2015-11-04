#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Create Sequence Coverage Plot.")

parser.add_argument("sect_profile_file", type=str,
                    help="The input profile file from KAT sect")

parser.add_argument("-p", "--output_type", type=str, default="png",
                    choices=["png", "pdf", "eps"],
                    help="The plot file type to create.")
parser.add_argument("-o", "--output", type=str, default="kat-profile",
                    help="The path to the output file.")
parser.add_argument("-t", "--title", type=str, default="Sequence Coverage Plot",
                    help="Title for plot")
parser.add_argument("-a", "--x_label", type=str, default="X",
                    help="Label for x-axis")
parser.add_argument("-b", "--y_label", type=str, default="Y",
                    help="Label for y-axis")
parser.add_argument("-x", "--x_max", type=int, default=1000,
                    help="Maximum value for x-axis")
parser.add_argument("-y", "--y_max", type=int, default=1000,
                    help="Maximum value for y-axis")
parser.add_argument("-w", "--width", type=int, default=1024,
                    help="Width of canvas")
parser.add_argument("-l", "--height", type=int, default=1024,
                    help="Height of canvas")
parser.add_argument("-n", "--index", type=int, default=0,
                    help="Index of fasta entry to plot")
parser.add_argument("-d", "--header", type=str,
                    help="Name of fasta entry to plot (has priority over index)")
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                    help="Print extra information")
parser.set_defaults(verbose=False)

args = parser.parse_args()
# ----- end command line parsing -----

