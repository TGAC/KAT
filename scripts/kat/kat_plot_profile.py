#!/usr/bin/env python3

import sys
import argparse
import math
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from kat_plot_misc import *

def main():

	# ----- command line parsing -----
	parser = argparse.ArgumentParser(
		description="Create Sequence Coverage Plot.")

	parser.add_argument("sect_profile_file", type=str,
						help="The input profile file from KAT sect")

	parser.add_argument("sect_profile_file_2", nargs="?", type=str,
						help="The optional second input profile file from KAT sect")

	parser.add_argument("-o", "--output", type=str, default="kat-profile",
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
	parser.add_argument("-c", "--y2_label", type=str,
						help="Label for second y-axis")
	parser.add_argument("-X", "--x_max", type=int,
						help="Maximum value for x-axis")
	parser.add_argument("-x", "--x_min", type=int,
						help="Minimum value for x-axis")
	parser.add_argument("-Y", "--y_max", type=int,
						help="Maximum value for y-axis")
	parser.add_argument("-y", "--y_min", type=int,
						help="Minimum value for y-axis")
	parser.add_argument("-z", "--y2_max", type=int,
						help="Maximum value for second y-axis")
	parser.add_argument("-w", "--width", type=int, default=8,
						help="Width of canvas")
	parser.add_argument("-l", "--height", type=int, default=2.5,
						help="Height of canvas")
	parser.add_argument("-n", "--index", type=str, default="0",
						help="Comma separate list of indexes of fasta entry " \
							 "to plot")
	parser.add_argument("-d", "--header", type=str,
						help="Name of fasta entry to plot (has priority " \
							 "over index)")
	parser.add_argument("--dpi", type=int, default=300,
						help="Resolution in dots per inch of output graphic.")
	parser.add_argument("-v", "--verbose", dest="verbose",
						action="store_true",
						help="Print extra information")
	parser.set_defaults(verbose=False)

	args = parser.parse_args()
	# ----- end command line parsing -----

	names = []
	profiles = {}
	names2 = []
	profiles2 = {}

	input_file = open(args.sect_profile_file)

	last_name = ""
	for line in input_file:
		if line[0] == '>':
			last_name = line[1:-1]
			names.append(last_name)
		else:
			profiles[last_name] = line[:-1]
	input_file.close()

	input2_file = None
	if args.sect_profile_file_2 is not None:
		input2_file = open(args.sect_profile_file_2)

		last_name = ""
		for line in input2_file:
			if line[0] == '>':
				last_name = line[1:-1]
				names2.append(last_name)
			else:
				profiles2[last_name] = line[:-1]
		input2_file.close()

	if args.sect_profile_file_2 is not None and len(names) != len(names2):
		print("First and second input files are not the same length", file=sys.stderr)
		exit(1)

	if args.header is not None:
		names = [args.header]
	else:
		indexes = list(map(int, args.index.split(',')))
		names = [names[i] for i in indexes]

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
		if args.sect_profile_file_2 is None:
			y_label = "Coverage"
		else:
			y_label = "Coverage (first file)"

	if args.y2_label is not None:
		y2_label = args.y2_label
	else:
		y2_label = "Coverage (second file)"

	fig, axs = plt.subplots(len(names), 1, figsize=(args.width, args.height * (len(names) + 0.3)))

	pstrs = [profiles[name] for name in names]
	profs = [np.fromstring(pstr, dtype=float, sep=' ') for pstr in pstrs]
	if args.x_max is not None:
		maxlen = args.x_max
	else:
		maxlen = max(list(map(len, profs)))

	if args.x_min is not None:
		minlen = args.x_min
	else:
		minlen = 1

	maxval1 = max(list(map(max, profs)))

	pstrs2 = []
	profs2 = []
	maxval2 = 0
	if args.sect_profile_file_2 is not None:
		pstrs2 = [profiles2[name] for name in names]
		profs2 = [np.fromstring(pstr2, dtype=float, sep=' ') for pstr2 in pstrs2]
		maxval2 = max(list(map(max, profs2)))

	for i in range(len(names)):
		if names[i] not in profiles:
			sys.exit("Entry {:s} not found.".format(names[i]))
		else:

			profile = profs[i]

			profile2 = None
			if args.sect_profile_file_2 is not None:
				profile2 = profs2[i]
				if len(profile) != len(profile2):
					print("First and second input files are not the same length", file=sys.stderr)
					exit(1)

			## axis fix
			if len(names) > 1:
				ax1 = axs[i]
			else:
				ax1 = axs

			ax2 = ax1.twinx()
			x = np.arange(1, len(profile) + 1)

			ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
			ax1.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
			ax1.set_xlim(minlen, maxlen + 1)

			if i == len(names) - 1:
				ax1.set_xlabel(x_label)
				for tick in ax1.get_xticklabels():
					tick.set_rotation(90)
					tick.set_visible(True)
			else:
				ax1.set_xlabel("")
				for tick in ax1.get_xticklabels():
					tick.set_rotation(90)
					tick.set_visible(False)

			## y limits from args or auto
			if args.y_max is not None:
				maxval1 = args.y_max
				maxval2 = args.y_max
			if args.y_min is not None:
				minval = args.y_min
			else:
				minval = 1

			ax1.set_title(names[i], fontsize=12)
			ax1.set_ylim(minval, maxval1 * 1.1)
			ax1.set_ylabel(y_label, color='r')
			ax1.plot(x, profile, 'r-')

			if profile2 is not None:
				ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
				ax2.set_ylim(minval, maxval2 * 1.1)
				ax2.set_ylabel(y2_label, color='b')
				ax2.plot(x, profile2, 'b-')

	plt.tight_layout()

	st = plt.suptitle(title, fontsize=18)
	st.set_y(0.95)
	plt.subplots_adjust(top=0.85)

	if args.output_type is not None:
		output_name = args.output + '.' + args.output_type
	else:
		output_name = args.output

	plt.savefig(correct_filename(output_name), dpi=args.dpi)

if __name__ == '__main__':
	main()
