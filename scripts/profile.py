#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Create Sequence Coverage Plot.")

parser.add_argument("sect_profile_file", type=str,
                    help="The input profile file from KAT sect")

parser.add_argument("-o", "--output", type=str, default="kat-profile",
                    help="The path to the output file.")
parser.add_argument("-p", "--output_type", type=str,
                    help="The plot file type to create (default is based on given output name).")
parser.add_argument("-t", "--title", type=str, default="Sequence Coverage Plot",
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
parser.add_argument("-l", "--height", type=int, default=3,
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

names = []
profiles = {}

input_file = open(args.sect_profile_file)

last_name = ""
for line in input_file:
    if line[0] == '>':
        last_name = line[1:-1]
        names.append(last_name)
    else:
        profiles[last_name] = line[:-1]

if args.header is not None:
    name = args.header
else:
    name = names[args.index]

if name not in profiles:
    sys.exit("Entry {:s} not found.".format(name))
else:
    pstr = profiles[name]

profile = np.fromstring(pstr, dtype=int, sep=' ')

x = np.arange(1,len(profile)+1)

plt.figure(num = None, figsize=(args.width, args.height))
plt.gcf().subplots_adjust(bottom=0.15) # makes room for x-label

plt.plot(x, profile)

plt.title(name)
plt.xlabel("Position")
plt.ylabel("Coverage")
plt.grid(True, color="black", alpha=0.2)

if args.output_type is not None:
    output_name = args.output + '.' + args.output_type
else:
    output_name = args.output

plt.savefig(output_name, dpi=300)
