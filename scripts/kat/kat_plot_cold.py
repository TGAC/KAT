#!/usr/bin/env python3

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import ScalarFormatter
import math
from scipy import stats

def main():

    # ----- command line parsing -----
    parser = argparse.ArgumentParser(
        description="Creates a scatter plot that shows points for each assembly contig, " \
            "that are sized by sequence length, coloured by assembly duplication level.  Each point is located " \
            "on a scatter plot with logscale read k-mer coverage on the Y-axis and contig GC% on the X.")

    parser.add_argument("stats_file", type=str,
                        help="The stats file produced by 'kat cold'")

    parser.add_argument("-o", "--output", type=str, default=None,
                        help="The path to the output file.")
    parser.add_argument("-p", "--output_type", type=str,
                        help="The plot file type to create (default is based on " \
                        "given output name).")
    parser.add_argument("-t", "--title", type=str,
                        help="Title for plot")
    parser.add_argument("-y", "--y_max", type=int,
                        help="Maximum value for y-axis")
    parser.add_argument("-w", "--width", type=int, default=8,
                        help="Width of canvas")
    parser.add_argument("-l", "--height", type=int, default=6,
                        help="Height of canvas")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Resolution in dots per inch of output graphic.")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true",
                        help="Print extra information")
    parser.set_defaults(verbose=False)

    args = parser.parse_args()
    # ----- end command line parsing -----

    if args.verbose:
        print("\nCoLD plotting:", args.stats_file)


    if args.title is not None:
        title = args.title
    else:
        title = "KAT Contig Length and Duplication plot"

    x_label = "GC%"
    y_label = "Median K-mer Coverage"

    sizes=[]
    gcs=[]
    dups=[]
    covs=[]
    with open(args.stats_file) as infile:
        for l in infile:
            line = l.strip()
            if not line or line=="" or line.startswith("seq_name"):
                continue

            parts = line.split("\t")
            sizes.append(int(parts[5]))
            gcs.append(float(parts[4]) * 100.0)
            dups.append(int(parts[3]))
            covs.append(float(parts[1]))

    #for i in range(len(sizes)):
    #    print(sizes[i], gcs[i], covs[i], dups[i])

    if args.verbose:
        print("Stats file loaded into memory.  Found", len(sizes), "contigs.")
        print()
        print("Median K-mer Coverage:", stats.describe(covs))
        print("GC%:", stats.describe(covs))
        print("Sequence lengthsstats:", stats.describe(sizes))
        print("Duplication levels:", stats.describe(dups))
        print()

    # Check dups
    for i, dup in enumerate(dups):
        if dup <= 0:
            raise ValueError("Found a duplication level of: " + dup + ".  We require duplications levels to be >= 1.")
        elif dup >= 7:
            dups[i] = 6


    colours = ["#ef292980",
               "#ad7fa880",
               "#8ae23480",
               "#729fcf80",
               "#f2c27e80",
               "#fcaf3e80",
               "#fce94f80"]

    # Make sure the ranges are at least 25 on both axis, less doesn't make much sense
    ymax = args.y_max if args.y_max else max(covs) * 5
    ymax = max(ymax, 25)

    if args.verbose:
        print("Axis limits:")
        print("ymax:", ymax)

    fig = plt.figure(figsize=(args.width, args.height))
    ax = fig.add_subplot(111)

    ax.set_xlim([0.0, 100.0])
    ax.set_ylim([0.9, float(ymax)])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    h = [None] * len(sizes)

    for i in range(len(sizes)):
        size = math.sqrt(sizes[i])
        h[i] = ax.scatter(gcs[i], covs[i], color=colours[dups[i]-1], marker="o", s=size, edgecolors='black')

    ax.xaxis.grid(True, which='major')
    ax.yaxis.grid(True, which='major')
    ax.set_axisbelow(True)

    ax.set_title(title)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ScalarFormatter())

    dupsleg = []
    for i in range(6):
        dupsleg.append(mpatches.Patch(color=colours[i], alpha=1))

    sizeleg = []
    for i in [1000, 10000, 100000, 1000000]:
        size = math.sqrt(math.sqrt(i))
        sizeleg.append(mlines.Line2D([0], [0], linestyle="none", marker="o", markersize=size, markeredgecolor="black", markerfacecolor="gray"))

    legend1 = ax.legend(dupsleg, ["1x", "2x", "3x", "4x", "5x", "6x+"], ncol=1, scatterpoints=1, fontsize="small", bbox_to_anchor=(1.15, 1.0))
    legend2 = ax.legend(sizeleg, ["1Kbp", "10Kbp", "100Kbp", "1Mbp"], ncol=4, markerscale=1, numpoints=1, scatterpoints=1, labelspacing=2,
                        handletextpad=1.5, borderaxespad=1.5, fontsize="small", loc="upper center")
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)

    if args.output:
        if args.output_type is not None:
            output_name = args.output + '.' + args.output_type
        else:
            output_name = args.output

        plt.savefig(output_name, dpi=args.dpi)
    else:
        plt.show()
    plt.close()

if __name__ == '__main__':
	main()
