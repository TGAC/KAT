#!/usr/bin/env python3

import argparse
import abc
import sys
import traceback
import time
import matplotlib.pyplot as plt

try:
	from spectra import KmerSpectra, GCSpectra
except:
	from kat.spectra import KmerSpectra, GCSpectra

version = "2.X.X"
try:
	import kat
	version = kat.__version__
except:
	pass



def plot_hist(h, xmax, ymax, label="", xlab='Kmer Frequency', ylab='# Distinct Kmers', color=None):
	plt.plot([min(ymax, x) for x in h[:xmax]], label=label, color=color)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.legend()



def plot_hist_df(self, h, xmax, ymax, xlab='Kmer Frequency', ylab='# Distinct Kmers'):
	plt.plot([max(-ymax, min(ymax, x)) for x in [h[i + 1] - h[i] for i in range(xmax)]])
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.legend()

class SpectraAnalysis(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, haploid=False, freq_cutoff=10000, hom_peak_freq=0, k=27):
		self.k = k
		self.haploid = haploid
		self.freq_cutoff = freq_cutoff
		self.hom_peak = hom_peak_freq
		self.limx = 0
		self.limy = 0

	@abc.abstractmethod
	def analyse(self, min_elements=1, verbose=False):
		pass

	@abc.abstractmethod
	def plot(self, points=0, cap=0, to_screen=False, to_files=None):
		pass

	@abc.abstractmethod
	def peak_stats(self):
		pass

class HistKmerSpectraAnalysis(SpectraAnalysis):
	def __init__(self, filename, haploid=False, freq_cutoff=10000, hom_peak_freq=0, k=27):
		SpectraAnalysis.__init__(self, haploid=haploid, freq_cutoff=freq_cutoff, hom_peak_freq=hom_peak_freq, k=k)
		self.spectra = KmerSpectra(self.read_hist(filename, freq_cutoff), haploid=haploid, k=k)


	def read_hist(self, name, freq_cutoff=10000):
		f = open(name)
		histogram = [int(x.split()[1]) for x in f.readlines() if x and x[0] != '#'][:freq_cutoff]
		f.close()
		return histogram

	def plot(self, xmax=0, ymax=0, to_screen=False, to_files=None):
		if 0 == xmax: xmax = self.limx
		if 0 == ymax: ymax = self.limy
		print()
		print("Creating plots")
		print("--------------")
		print()

		if len(self.spectra.peaks) == 0:
			print("No peaks in K-mer frequency histogram.  Not plotting.")
		else:
			print("Plotting K-mer frequency distributions...")
			#self.spectra.printPeaks()

			fig = plt.figure()
			plot_hist(self.spectra.histogram, xmax, ymax, label="Actual", color="black")

			for p_i, p in enumerate(self.spectra.peaks, start=1):
				plot_hist(p.Ty, xmax, ymax, label=p.description)

			plot_hist(self.spectra.Ty, xmax, ymax, label="Fitted model", color="green")

			if to_screen:
				plt.show()

			if to_files:
				filename = to_files + ".dists.png"
				fig.savefig(filename)
				print("- Saved plot to:", filename)

		print()

	def analyse(self, min_elements=1, verbose=False):

		# Create the peaks
		if verbose:
			print("Analysing spectra")
		self.spectra.analyse(min_elements=min_elements, verbose=verbose)
		if self.spectra.peaks:
			self.limy = int(max(int(self.spectra.maxValue() * 1.1 / 1000) * 1000, self.limy))
			self.limx = int(max(min(self.spectra.peaks[-1].mean() * 2, len(self.spectra.histogram)), self.limx))

		if verbose:
			print("Plot limits: y->%d, x->%d" % (self.limy, self.limx))

	def peak_stats(self):
		print()
		print("K-mer frequency spectra statistics")
		print("----------------------------------")
		self.spectra.printGenomeStats(self.hom_peak)

class GCKmerSpectraAnalysis(SpectraAnalysis):
	def __init__(self, filename, haploid=False, freq_cutoff=10000, hom_peak_freq=0, k=27):
		SpectraAnalysis.__init__(self, haploid=haploid, freq_cutoff=freq_cutoff, hom_peak_freq=hom_peak_freq, k=k)
		cov_histo, gc_histo = self.read_file(filename, freq_cutoff)
		self.mean_gc = sum([i * x for i, x in enumerate(gc_histo)]) / sum(gc_histo)
		self.cov_spectra = KmerSpectra(cov_histo, haploid=haploid, k=k)
		self.gc_dist = GCSpectra(gc_histo, k=k)		# Not really a Kmer spectra but let's shoehorn it in anyway


	def read_file(self, name, freq_cutoff=10000, k=27):
		f = open(name)
		cov_histogram = None
		gc_histogram = []
		for x in f.readlines():
			if x and x[0] != '#':
				# Kmer coverage histo
				parts = x.split()
				gc_histogram.append(sum([int(y) for y in parts]))
				if not cov_histogram:
					cov_histogram = [0] * len(parts)
				for i, y in enumerate(parts):
					cov_histogram[i] += int(y)
		f.close()
		return cov_histogram, gc_histogram


	def plot(self, xmax=0, ymax=0, to_screen=False, to_files=None):
		if 0 == xmax: xmax = self.limx
		if 0 == ymax: ymax = self.limy

		print()
		print("Creating plots")
		print("--------------")
		print()
		if len(self.cov_spectra.peaks) == 0:
			print("No peaks in K-mer frequency histogram.  Not plotting.")
		else:

			print("Plotting K-mer frequency distributions...")

			fig = plt.figure()
			plot_hist(self.cov_spectra.histogram, xmax, ymax, label="Actual", color="black")

			for p_i, p in enumerate(self.cov_spectra.peaks, start=1):
				plot_hist(p.Ty, xmax, ymax, label=p.description)

			plot_hist(self.cov_spectra.Ty, xmax, ymax, label="Fitted model", color="green")

			if to_screen:
				plt.show()

			if to_files:
				filename = to_files + ".dists.png"
				fig.savefig(filename)
				print(" - Saved plot to:", filename)

		if len(self.gc_dist.peaks) == 0:
			print("No peaks in GC distribution.  Not plotting.")
		else:
			print("Plotting GC distributions...")
			#self.gc_dist.printPeaks()

			xmax = self.gc_dist.k
			ymax = max(self.gc_dist.histogram) * 1.1

			fig = plt.figure()
			plot_hist(self.gc_dist.histogram, xmax, ymax, label="Actual", xlab="GC count", color="black")

			for p_i, p in enumerate(self.gc_dist.peaks, start=1):
				plot_hist(p.Ty, xmax, ymax, label="Dist %d" % p_i, xlab="GC count")

			plot_hist(self.gc_dist.Ty, xmax, ymax, label="Fitted model", xlab="GC count", color="green")

			if to_screen:
				plt.show()

			if to_files:
				filename = to_files + ".gc.png"
				fig.savefig(filename)
				print(" - Saved plot to:", filename)

		print()


	def analyse(self, min_elements=1, verbose=False):

		# Create the peaks
		if verbose:
			print("Analysing K-mer spectra")
		self.cov_spectra.analyse(min_elements=min_elements, verbose=verbose)
		if self.cov_spectra.peaks:
			self.limy = int(max(int(self.cov_spectra.maxValue() * 1.1 / 1000) * 1000, self.limy))
			self.limx = int(max(min(self.cov_spectra.peaks[-1].right() * 1.1, len(self.cov_spectra.histogram)), self.limx))

		if verbose:
			print()
			print("Plot limits: y->%d, x->%d" % (self.limy, self.limx))
			print()

		# Create the peaks
		if verbose:
			print("Analysing GC distribution")
		self.gc_dist.analyse(min_elements=min_elements, verbose=verbose)


	def peak_stats(self):
		print()
		print("K-mer frequency spectra statistics")
		print("----------------------------------")
		self.cov_spectra.printGenomeStats(self.hom_peak)
		print()
		print("GC distribution statistics")
		print("--------------------------")
		# Step 4, genome stats
		print("K-value used:", str(self.gc_dist.k))
		print("Peaks in analysis:", str(len(self.gc_dist.peaks)))
		print()
		self.gc_dist.printPeaks()
		print("Mean GC:", "{0:.2f}".format(self.mean_gc * (100.0 / self.gc_dist.k)) + "%")




class MXKmerSpectraAnalysis(SpectraAnalysis):
	def __init__(self, filename, cns_cutoff=3, haploid=False, freq_cutoff=10000, hom_peak_freq=0, k=27):
		SpectraAnalysis.__init__(self, haploid=haploid, freq_cutoff=freq_cutoff, hom_peak_freq=hom_peak_freq, k=k)
		self.spectras = [KmerSpectra(self.read_mx(filename, freq_cutoff=freq_cutoff, column=0, cumulative=True), haploid=haploid, k=k)]
		for i in range(cns_cutoff):
			self.spectras.append(KmerSpectra(self.read_mx(filename, freq_cutoff=freq_cutoff, column=i, cumulative=False), haploid=haploid, k=k))

	def read_mx(self, name, freq_cutoff=10000, column=1, cumulative=False):
		f = open(name)

		histogram = []
		if cumulative:
			histogram = [sum([int(y) for y in x.split()[column:]]) for x in f.readlines() if x and x[0] != '#'][
							 :freq_cutoff][1:]
		else:
			histogram = [int(x.split()[column]) for x in f.readlines() if x and x[0] != '#'][:freq_cutoff][1:]
		f.close()
		return histogram

	def plot(self, xmax=0, ymax=0, to_screen=False, to_files=None):
		if 0 == xmax: xmax = self.limx
		if 0 == ymax: ymax = self.limy
		print()
		print("Creating plots")
		print("--------------")
		print()
		for s_i, s in enumerate(self.spectras):
			if self.spectras[0] == s:
				slabel = "General Spectra"
			else:
				slabel = "%dx present" % (s_i - 1)
			fig = plt.figure()

			print("Plotting:",slabel)
			s.printPeaks()

			plot_hist(s.histogram, xmax, ymax, label=slabel, color="black")
			plot_hist(s.Ty, xmax, ymax, label=slabel + " fit")
			for p_i, p in enumerate(s.peaks, start=1):
				plot_hist(p.Ty, xmax, ymax, label=p.description)

			if to_screen:
				plt.show()

			if to_files:
				suffix = ".spectra-cn" + str(s_i-1) + ".png"
				if self.spectras[0] == s:
					suffix = ".general.png"
				filename = to_files + suffix
				fig.savefig(filename)
				print(" - Saved plot to:", filename)

			print()


	def analyse(self, min_elements=1, verbose=False):

		maxValue = 0
		right = 0
		for s_i, s in enumerate(self.spectras):
			if s_i == 0:
				print("\nAnalysing full spectra")
			else:
				print("\nAnalysing spectra with copy number", s_i-1)
			s.analyse(min_elements=min_elements, verbose=verbose)
			if s.peaks:
				maxValue = max(maxValue, s.maxValue())
				right = max(right, s.peaks[-1].right())
			elif s_i == 0:
				raise RuntimeError("No peaks detected for full spectra.  Can't continue.")

		self.limy = int(max(int(maxValue * 1.1 / 1000) * 1000, self.limy))
		self.limx = int(max(min(right * 1.1, len(s.histogram)), self.limx))

		print("\nAnalysed spectra for all requested copy numbers.")

		if verbose:
			print("\nPlot limits: y->%d, x->%d" % (self.limy, self.limx))

	def peak_stats(self):
		"""TODO: Runs analyse (TODO:include general spectra)
				 Takes enough peaks as to cover a given % of the elements:
					 - Find the peak across all distributions
					 - Reports peak stats
				 If multiple peaks have been analyzed, tries to find the "main unique" and explains the results based on that freq.
		"""
		# First check to see if we have anything to work with
		if len(self.spectras[0].peaks) == 0:
			raise ValueError("Main spectra distribution does not contain any peaks.")

		# step 1, try to find a reasonable mean for kmer frequency.
		# weighted means by number of elements?
		print("\nMain spectra statistics")
		print("-----------------------")
		self.spectras[0].printGenomeStats(self.hom_peak)

		# step 2, try to estimate the assembly completeness
		completeness = self.calcAssemblyCompleteness()
		print("Estimated assembly completeness:", ("{0:.2f}".format(completeness) + "%") if completeness > 0.0 else "Unknown")

		# step 3, selects frequencies for peaks from bigger to smaller till X% of the elements are covered or no more peaks
		print("\nBreakdown of copy number composition for each peak")
		print("--------------------------------------------------")

		general_dists = self.spectras[0].peaks
		goal = 0.99 * sum([x.elements() for x in general_dists])
		maxpeaks = 10
		general_dists.sort(key=lambda x: -x.elements())
		af = []
		peaks = 0
		covered = 0
		for x in general_dists:
			af.append(x.mean())
			peaks += 1
			covered += x.elements()
			if peaks == maxpeaks or covered > goal:
				break

		# step 3, report for each peak
		# get the candidate peak on each spectra
		for f in af:
			total = 0
			pd = {}
			for i, s in enumerate(self.spectras, start=1):
				m = [(x.mean(), x.elements()) for x in s.peaks if x.mean() > 0.8 * f and x.mean() < 1.2 * f]
				if len(m) == 1:
					pd[i] = m[0]
					total += m[0][1]
				if len(m) > 1:
					print("WARNING, MORE THAT 1 PEAK FOR f=%.3f FOUND ON THE %dx SPECTRA!!!" % (f, i))
			print("\n---- Report for f=%.3f (total elements %d)----" % (f, total))
			for i in range(len(self.spectras) - 1):
				if i in list(pd.keys()):
					print(
						" %dx: %.2f%% (%d elements at f=%.2f)" % (i, float(pd[i][1]) * 100 / total, pd[i][1], pd[i][0]))
				else:
					print(" %dx: No significant content" % i)


	def calcAssemblyCompleteness(self):

		if self.spectras[0].peaks:
			hpi = self.spectras[0].getHomozygousPeakIndex(self.hom_peak)
			opt_freq = int(self.spectras[0].peaks[hpi-1].mean())
			absent_count = self.spectras[1].histogram[opt_freq]
			present_count = self.spectras[2].histogram[opt_freq]
			return (present_count / (absent_count + present_count)) * 100.0
		else:
			return 0.0

def get_properties_from_file(input_file):
	k = 27
	mx = False
	gcp = False

	f = open(input_file)
	i = 0
	for l in f.readlines():
		if i > 10:
			break
		line = l.strip()
		if line.startswith("#"):
			if line.startswith("# Kmer value:"):
				k = int(line.split(":")[1])
			elif line.startswith("# Rows:"):
				mx = True
			elif line.startswith("# YLabel:GC count"):
				gcp = True
		i+=1
	f.close()

	return k, mx, gcp

def main():
	# ----- command line parsing -----
	parser = argparse.ArgumentParser(
		description="""Analyse a comp matrix file with respect to the distributions and copy numbers seen within.""")

	parser.add_argument("input", type=str,
						help="The input should be either a KAT spectra-cn matrix file a KAT GCP matrix file or a KAT histogram file.")

	parser.add_argument("-o", "--output_prefix",
						help="If present then plots are sent to files starting with this prefix.")
	parser.add_argument("-c", "--cns", type=int, default=4,
						help="The number of copy numbers to consider in the analysis.  Only applicable if input is a spectra-cn matrix file.")
	parser.add_argument("-f", "--freq_cutoff", type=int, default=500,
						help="The maximum frequency cutoff point to consider.  Analysis will be done up to this frequency.")
	parser.add_argument("-e", "--min_elem", type=int, default=10000,
						help="Any new distribution that adds less to this number of distinct K-mers will not be added.")
	parser.add_argument("-p", "--plot", action='store_true',
						help="Plot best cumulative fit for all peaks.")
	parser.add_argument("-z", "--homozygous_peak", type=int, default=0,
						help="The approximate kmer frequency for the homozygous peak.  Allows us to calculate a more accurate genome size estimate.")
	parser.add_argument("--haploid", action='store_true',
						help="If selected then we do not try to detect a heterozygous peak")
	parser.add_argument("-v", "--verbose", action='store_true',
						help="Print additional information.")
	parser.add_argument("--from_kat", action='store_true', help=argparse.SUPPRESS)

	args = parser.parse_args()

	if not args.from_kat:
		print("KAT K-mer Distribution Analysis Script")
		print("Version:", version)
		print()
		if args.verbose:
			print("Analysing distributions for:", args.input)
		else:
			print("Analysing distributions for:", args.input, "... ", end="", flush=True)

	k, mx, gcp = get_properties_from_file(args.input)
	if args.verbose:
		print("Input file generated using K", k)

	a = None
	if mx:
		if gcp:
			if args.verbose:
				print("GC vs Coverage matrix file detected")
			a = GCKmerSpectraAnalysis(args.input, haploid=args.haploid, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)
		else:
			if args.verbose:
				print("Copy number spectra matrix file detected")
				print("Processing", args.cns, "spectra")
			a = MXKmerSpectraAnalysis(args.input, haploid=args.haploid, cns_cutoff=args.cns, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)
	else:
		if args.verbose:
			print("Kmer coverage histogram file detected")
		a = HistKmerSpectraAnalysis(args.input, haploid=args.haploid, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)

	if not a:
		raise RuntimeError("Couldn't generate a valid spectra analysis")

	try:
		start = time.time()
		a.analyse(min_elements=args.min_elem, verbose=args.verbose)
		end = time.time()
		print(("" if args.verbose else "done.  ") + "Time taken: ", '{0:.1f}'.format(end - start) + 's')
		a.peak_stats()
		if args.plot or args.output_prefix:
			a.plot(xmax=args.freq_cutoff, to_screen=args.plot, to_files=args.output_prefix)
	except Exception:
		print("\nERROR\n-----", file=sys.stderr)
		traceback.print_exc(file=sys.stderr)

if __name__ == '__main__':
	main()
