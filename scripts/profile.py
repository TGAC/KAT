#!/usr/bin/env python

import sys
import argparse
import math
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
parser.add_argument("-l", "--height", type=int, default=3,
                    help="Height of canvas")
parser.add_argument("-n", "--index", type=str, default="0",
                    help="Comma separate list of indexes of fasta entry to plot")
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
input_file.close()

if args.header is not None:
    names = [args.header]
else:
    indexes = map(int, args.index.split(','))
    names = map(lambda i: names[i], indexes)

if args.title is not None:
    title = args.title
else:
    title = "Sequence Coverage Plot"

if args.x_label is not None:
    x_label = args.x_label
else:
    x_label = "Position"

if args.y_label is not None:
    y_label = args.y_label
else:
    y_label = "Coverage"

plt.figure(1, figsize=(args.width, args.height * len(names)))

pstrs = map(lambda name: profiles[name], names)
profs = map(lambda pstr: np.fromstring(pstr, dtype=int, sep=' '), pstrs)
maxlen = max(map(len, profs))
tickdist = maxlen/7
tickdist = int(round(tickdist, -int(math.floor(math.log10(tickdist)))))
for i in range(len(names)):
    if names[i] not in profiles:
        sys.exit("Entry {:s} not found.".format(names[i]))
    else:
        # pstr = profiles[names[i]]

        # profile = np.fromstring(pstr, dtype=int, sep=' ')
        profile = profs[i]

        x = np.arange(1,len(profile)+1)

        plt.subplot(len(names), 1, i+1)
        plt.plot(x, profile)
        plt.ylim(0,np.max(profile)*1.1)
        # diff = maxlen - len(profile)
        # plt.xlim(1 - diff/2, len(profile)+diff/2+1)
        plt.xlim(1, maxlen+1)
        xticks = range(0, len(profile)+1, tickdist)
        xticks[0] = 1
        if len(profile) - xticks[-1] > tickdist/2:
            xticks.append(len(profile))
        else:
            xticks[-1] = len(profile)
        plt.xticks(xticks)

        plt.title(names[i], fontsize=12)
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.grid(True, color="black", alpha=0.2)
        plt.tight_layout()

if title != "":
    plt.suptitle(title, fontsize=14)
    plt.subplots_adjust(top=0.95)

if args.output_type is not None:
    output_name = args.output + '.' + args.output_type
else:
    output_name = args.output

plt.savefig(output_name, dpi=300)
