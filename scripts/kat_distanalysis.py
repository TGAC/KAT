#!/usr/bin/env python3

import argparse
import abc
import sys

import numpy as np
from scipy import mean, optimize
import matplotlib.pyplot as plt


def plot_hist(h, points, cap, label=""):
	plt.plot([min(cap, x) for x in h[:points]], label=label)


class KmerSpectra(object):
	"""A kmer spectra, comprised of different peaks.
	Contains the general fitting method"""

	###----------------- INIT, LOAD, DUMP, ETC --------------------
	def __init__(self, histogram, k=27):
		"""Init with the objective histogram"""
		self.histogram = histogram
		self.peaks = []
		self.k = k

	def total_values(self, start=1, end=10000):
		return list(map(sum, list(zip(*[x.points(start, end) for x in self.peaks]))))

	###----------------- FIND AND CREATE PEAKS --------------------


	def analyse(self,min_perc=1, min_elem=100000, verbose=False):
		if verbose:
			print("Creating initial peaks...", end="")
			sys.stdout.flush()
		self.create_peaks(min_perc=min_perc, min_elem=min_elem, verbose=verbose)

		if self.peaks:
			if verbose:
				print("done.", len(self.peaks), "peaks initially created. fmin=%d (%d) | fmax=%d (%d)" % (self.fmin,self.histogram[self.fmin],self.fmax,self.histogram[self.fmax]))
			self.optimize_peaks(min_perc=min_perc, min_elem=min_elem, verbose=verbose)
			if verbose:
				print("Optimising all peaks in spectra...", end="")
				sys.stdout.flush()
			self.optimize_overall()
			if verbose:
				print("done.", len(self.peaks), "peaks present after optimisation:")
				self.printPeaks()
		elif verbose:
			print("done. No peaks created")

	def deriv(self, histogram):
		return [histogram[i + 2] - histogram[i] for i in range(len(histogram) - 2)]

	def progsmoothderiv(self, histo):
		''' Smoothed df/dx (dx=2 units)'''
		# always return mean
		smoothed = []
		for i in range(len(histo) - 2):
			delta = max(int(i/10), 1) if i > 0 else 1

			r1_start = i + 1
			r1_end = i + 1 + delta
			m1 = mean(histo[r1_start:r1_end])

			r2_start = max(i - delta, 0)
			r2_end = i+1
			m2 = mean(histo[r2_start:r2_end])

			smoothed.append(m1 - m2)
		return smoothed

	def find_maxima(self, center, radius, min_perc, min_elem, histo=None):

		start=center-radius
		end=center+radius

		hs = histo[start:end] if histo else self.histogram[start:end]

		if len(hs) == 0:
			raise ValueError("Histogram has length of 0: Center:", center, "; Radius:", radius)

		# Smoothed df/dx (dx=2 units)
		deriv = self.progsmoothderiv(hs)
		fmax = hs.index(max(hs))
		if fmax == 0 or fmax == len(hs) - 1: return 0

		# find a single inflection point-> that is the maxima
		# Reject by voting to avoid oscillations.
		failpoints = 0
		for i in range(fmax):
			if deriv[i] < 0: failpoints += 1
		for i in range(fmax + 1, len(deriv)):
			if deriv[i] > 0: failpoints += 1
		if float(failpoints) / (2 * radius + 1) > .1:
			return 0

		# TODO: discard maxima if not 1% of the already contained elements on peaks are there
		if sum(hs) < min((float(min_perc) / 100) * sum([x.elements for x in self.peaks]), min_elem):
			print("Distribution on %d too small to be considered (%d elements)" % (center - radius + fmax, sum(hs)))
			return 0
		# TODO: do further validation
		# print "maxima found on %d, for %s" % (fmax,hs)
		return center - radius + fmax

	def add_peak_and_update_cuts(self, lm, reset_opt=False):
		# Receives a local maxima, adds a peak there if there was none, updates the cuts.
		# If reset_opt is set, it will reset all the peak counts/sd and only keep the means
		fdists = [x.mean for x in self.peaks]
		for f in fdists:
			if lm >= f - f / 5 and lm <= f + f / 5:
				print("WARNING!!! Dist on %d is not to be added, because existing dist on %d is too close to it." % (
				lm, f))
				return False
		fdists.append(lm)
		fdists.sort()
		self.cuts = [self.fmin]
		for i in range(len(fdists) - 1):
			self.cuts.append(int(fdists[i] * (1 + float(fdists[i]) / (fdists[i] + fdists[i + 1]))))
		self.cuts.append(int(min(len(self.histogram) - 1, fdists[-1] * 1.5)))
		# print "cuts on %s" % str(self.cuts)
		if reset_opt:
			self.peaks = []
		for i in range(len(fdists)):
			if reset_opt or i == fdists.index(lm):
				self.peaks.append(KmerPeak(lm, fdists[i], fdists[i] / 6,
										   sum([self.histogram[j] for j in range(self.cuts[i], self.cuts[i + 1])]),
										   self.histogram[fdists[i]]))
				self.peaks.sort(key=lambda x: x.mean)
		return True

	def create_peaks(self, min_perc, min_elem, verbose=False):
		fmin = 0
		# walk till first local minimum (d(f)>0)
		for i in range(len(self.histogram)):
			if self.histogram[i] < self.histogram[i + 1]:
				fmin = i
				break
		if not fmin:
			print()
			print(self.histogram[:25])
			print("Could not find local minima")
			self.fmax = 0
			self.maxval = 0

		# Find the global fm|p(fm)=max(p)
		fmax = self.histogram.index(max(self.histogram[fmin:]))

		if fmax < 10:
			print("fmax<10, no analysis can be performed")
			return
		self.fmax = fmax
		self.fmin = fmin
		self.maxval = self.histogram[fmax]
		# Explore fm/x and fm*x for x in [1,4]
		for f in [fmax / 4.0, fmax / 3.0, fmax / 2.0, fmax, fmax * 2, fmax * 3, fmax * 4]:
			if int(f / 5.0) > 0 and f + f / 5.0 < len(self.histogram):
				lm = self.find_maxima(int(f), int(f / 5.0), min_perc, min_elem)
				if lm:
					self.add_peak_and_update_cuts(lm, reset_opt=True)
					# if not fdists:
					#    print "Local maxima not found, skipping further analysis"
					#    return
					# Guess counts for all relevant distributions
					# print "Local maxima found on %s" %(str(fdists))

	###----------------- OPTIMIZATION --------------------

	def optimize_peaks(self, min_perc, min_elem, verbose=False):
		sortedpeaks = [x for x in self.peaks]
		sortedpeaks.sort(key=lambda x: -x.elements)

		# create a base as the histogram and start from there
		base = [x for x in self.histogram]
		for p_i, p in enumerate(sortedpeaks):
			i = self.peaks.index(p)
			# locally optimize the preak
			if verbose:
				print("Optimizing peak:", p_i)
				print(" - current:", p, "at [%d,%d]" % (self.cuts[i], self.cuts[i + 1]))
			p.constrained_opt(self.cuts[i], base[self.cuts[i]:self.cuts[i + 1]])
			if verbose:
				print(" - new    :", p)
			# substract the peak effect from the baseline
			for i in range(len(self.histogram)):
				base[i] = base[i] - p.point(i)
		updated = False
		for f in [self.fmax / 4, self.fmax / 3, self.fmax / 2, self.fmax, self.fmax * 2, self.fmax * 3, self.fmax * 4]:
			if int(f + f / 5) < len(self.histogram) and sum(
					[base[x] for x in range(int(f - f / 5), int(f + f / 5))]) > .01 * sum(
					[p.elements for p in self.peaks]):
				lm = self.find_maxima(int(f), int(f / 5), min_perc, min_elem, base)
				if lm:
					if self.add_peak_and_update_cuts(lm, reset_opt=False):
						updated = True
						# else:
						#    print "No new maxima found on %d +/- %d" % (f,f/5)
		if updated: self.optimize_peaks(min_perc, min_elem)

	def residues(self, p):
		"""p has a triplet for each distribution"""
		# print p
		for i in range(len(self.peaks)):
			self.peaks[i].mean = p[i * 3]
			self.peaks[i].shape = p[i * 3 + 1]
			self.peaks[i].elements = p[i * 3 + 2]
		tv = self.total_values()
		r = [tv[i] - self.histogram[i] for i in range(self.cuts[0], self.cuts[-1])]
		# a lot more penalty for being UP:
		for i in range(len(r)):
			if r[i] > 0:
				r[i] = 2 * r[i]
		# for j in xrange(len(self.peaks)):
		#    if p[0]>self.peaks[j].target_max*1.1 or p[0]<self.peaks[j].target_max*.9:
		#        #for i in xrange(len(r)): r[i]+=r[i]*10*(float(p[j*3]-self.peaks[j].target_max)/self.peaks[j].target_max)
		#        for i in xrange(len(r)): r[i]+=r[i]*2
		return r

	def optimize_overall(self):
		p = []
		for pk in self.peaks:
			p.append(pk.mean)
			p.append(pk.shape)
			p.append(pk.elements)
		optimize.leastsq(self.residues, p)
		# once the better fit is found, check if by taking the unfitted elements new distributions arise.
		return


	def getHomozygousPeakIndex(self, approx_freq=0):
		# Work out which peak is the homozygous peak
		peak_index = len(self.peaks)
		if approx_freq > 0:
			min_peak_index = peak_index
			delta_peak_freq = 1000000
			for p_i, p in enumerate(self.peaks, start=1):
				delta = abs(p.mean - approx_freq)
				if delta_peak_freq > delta:
					delta_peak_freq = delta
					min_peak_index = p_i
			peak_index = min_peak_index

		return peak_index


	def calcGenomeSize(self, hom_peak=0):

		hom_peak_index = len(self.peaks) if hom_peak == 0 else hom_peak

		# Just use first spectra for this calculation
		sum = 0
		for p_i, p in enumerate(self.peaks, start=1):
			sum += p_i * p.elements

		return int(sum / hom_peak_index)


	def calcHetRate(self, genome_size=0, hom_peak=0):

		genomesize = genome_size if genome_size > 0 else self.calcGenomeSize()

		hom_peak_index = len(self.peaks) if hom_peak == 0 else hom_peak

		sum = 0
		if hom_peak_index < 2:
			return 0.0

		for p_i, p in enumerate(self.peaks, start=1):
			#Skip the last peak
			if p_i >= hom_peak_index:
				break

			sum += p.elements / self.k

		return (sum / genomesize) * 100.0


	def calcKmerCoverage(self):
		tot_vol = sum([x.elements for x in self.peaks])
		weighted = sum([x.mean * x.elements for x in self.peaks])
		return int(weighted / tot_vol)

	def printPeaks(self):
		if len(self.peaks) > 0:
			print("Index\t" + KmerPeak.getTabHeader())
			for p_i, p in enumerate(self.peaks, start=1):
				print(str(p_i) + "\t" + p.toTabString())
		else:
			print("No peaks detected")

	def printGenomeStats(self, hom_peak_freq):
		# Step 4, genome stats
		print("K-value used:", self.k)
		print("Peaks in analysis:", len(self.peaks))
		self.printPeaks()

		print("Mean k-mer frequency:", str(self.calcKmerCoverage()) + "x")

		hp = self.getHomozygousPeakIndex(hom_peak_freq)
		print("Homozygous peak index:", hp)
		gs = self.calcGenomeSize(hom_peak=hp)
		print("Estimated genome size:", '{0:.2f}'.format(float(gs) / 1000000.0), "Mbp")
		if (hp > 1):
			hr = self.calcHetRate(gs)
			print("Estimated heterozygous rate:", '{0:.2f}'.format(hr), "%")
		

class KmerPeak(object):
	"""A distribution representing kmers covered a certain number of times.
	Contains methods for fitting to an interval"""

	def __init__(self, target_max, mean, shape, elements, peak):
		self.target_max = target_max
		self.mean = mean
		self.shape = shape

		self.elements = elements
		self.peak = peak

	def __str__(self):
		return "<Kmer distribution peak of %d at %dx, with volume of %d elements>" % (self.peak, self.mean, self.elements)

	def toTabString(self):
		return "\t".join([str(int(self.mean)), str(int(self.peak)), str(int(self.elements))])

	@staticmethod
	def getTabHeader():
		return "Freq\tMax\tVolume"

	def point(self, x):
		"""Normalized Gaussian"""
		return float(self.elements) / np.sqrt(2 * np.pi) / self.shape * np.exp(
			-(x - self.mean) ** 2 / 2. / self.shape ** 2)

	def points(self, start, end):
		return [self.point(x) for x in range(start, end)]

	def residues(self, p, offset, histogram):
		"""residues of fitting self with parameters p to offset obj[0] with values obj[1]"""
		# TODO: enforce a constrain?
		self.mean = p[0]
		self.shape = p[1]
		self.elements = p[2]
		r = [self.point(offset + i) - histogram[i] for i in range(len(histogram))]
		# a lot more penalty for being UP:
		for i in range(len(r)):
			if r[i] > 0:
				r[i] = 2 * r[i]
				# if p[0]>self.target_max*1.1 or p[0]<self.target_max*.9:
				# for i in xrange(len(r)): r[i]+=r[i]*10*(float(p[0]-self.target_max)/self.target_max)
		# for i in xrange(len(r)): r[i]+=r[i]*2
		# penalizatin for moving the mean

		return r

	def constrained_opt(self, offset, histogram):
		"""Does a constrained optimization of fitting to match the histogram"""
		# print "Optimizing peak %s" % str(self)
		optimize.leastsq(self.residues, (self.mean, self.shape, self.elements), (offset, histogram))
		# self.mean=p[0]
		# self.shape=p[1]
		# self.elements=p[2]
		pass

def plot_hist(h, points, cap, label=""):
	plt.plot([min(cap, x) for x in h[:points]], label=label)
	plt.xlabel('Kmer Frequency')
	plt.ylabel('# Distinct Kmers')
	plt.legend()

def plot_hist_df(self, h, points, cap):
	plt.plot([max(-cap, min(cap, x)) for x in [h[i + 1] - h[i] for i in range(points)]])
	plt.xlabel('Kmer Frequency')
	plt.ylabel('# Distinct Kmers')
	plt.legend()

class SpectraAnalysis(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, freq_cutoff=10000, hom_peak_freq=0, k=27):
		self.k = k
		self.freq_cutoff = freq_cutoff
		self.hom_peak = hom_peak_freq
		self.limx = 0
		self.limy = 0

	@abc.abstractmethod
	def analyse(self, min_perc=1, min_elem=100000, verbose=False):
		pass

	@abc.abstractmethod
	def plot(self, points=0, cap=0):
		pass

	@abc.abstractmethod
	def peak_stats(self):
		pass

class HistKmerSpectraAnalysis(SpectraAnalysis):
	def __init__(self, filename, freq_cutoff=10000, hom_peak_freq=0, k=27):
		SpectraAnalysis.__init__(self, freq_cutoff=freq_cutoff, hom_peak_freq=hom_peak_freq, k=k)
		self.spectra = KmerSpectra(self.read_hist(filename, freq_cutoff), k=k)


	def read_hist(self, name, freq_cutoff=10000):
		f = open(name)
		histogram = [int(x.split()[1]) for x in f.readlines() if x and x[0] != '#'][:freq_cutoff]
		f.close()
		return histogram

	def plot(self, points=0, cap=0):
		if 0 == points: points = self.limx
		if 0 == cap: cap = self.limy
		print()
		print("Creating plot")
		print("-------------")
		print()
		self.spectra.printPeaks()
		print()

		plt.figure()
		plot_hist(self.spectra.histogram, points, cap, label="Histogram")
		plot_hist(self.spectra.total_values(1, points + 1), points, cap, label="Fitted distribution")

		for p_i, p in enumerate(self.spectra.peaks, start=1):
			plot_hist(p.points(1, points + 1), points, cap, label="fit dist %d" % p_i)
		plt.show()

	def analyse(self, min_perc=1, min_elem=100000, verbose=False):

		# Create the peaks
		print("\nAnalysing spectra")
		self.spectra.analyse(min_perc=min_perc, min_elem=min_elem, verbose=verbose)
		if self.spectra.peaks:
			self.limy = int(max(int(self.spectra.maxval * 1.1 / 1000) * 1000, self.limy))
			self.limx = int(max(min(self.spectra.peaks[-1].mean * 2, len(self.spectra.histogram)), self.limx))

		if verbose:
			print("\nPlot limits: y->%d, x->%d" % (self.limy, self.limx))

	def peak_stats(self):
		print()
		print("Spectra statistics")
		print("------------------")
		self.spectra.printGenomeStats(self.hom_peak)

class MXKmerSpectraAnalysis(SpectraAnalysis):
	def __init__(self, filename, cns_cutoff=3, freq_cutoff=10000, hom_peak_freq=0, k=27):
		SpectraAnalysis.__init__(self, freq_cutoff=freq_cutoff, hom_peak_freq=hom_peak_freq, k=k)
		self.spectras = [KmerSpectra(self.read_mx(filename, freq_cutoff=freq_cutoff, column=0, cumulative=True), k=k)]
		for i in range(cns_cutoff):
			self.spectras.append(KmerSpectra(self.read_mx(filename, freq_cutoff=freq_cutoff, column=i, cumulative=False), k=k))

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

	def plot(self, points=0, cap=0):
		if 0 == points: points = self.limx
		if 0 == cap: cap = self.limy
		print()
		print("Creating plots")
		print("--------------")
		print()
		for s_i, s in enumerate(self.spectras):
			if self.spectras[0] == s:
				slabel = "General Spectra"
			else:
				slabel = "%dx present" % (s_i - 1)
			plt.figure()

			print("Plotting:",slabel)
			s.printPeaks()
			print()

			plot_hist(s.histogram, points, cap, label=slabel)
			plot_hist(s.total_values(1, points + 1), points, cap, label=slabel + " fit")
			for p_i, p in enumerate(s.peaks, start=1):
				plot_hist(p.points(1, points + 1), points, cap, label="fit dist %d" % p_i)

			plt.show()

	def analyse(self, min_perc=1, min_elem=100000, verbose=False):

		for s_i, s in enumerate(self.spectras):
			if s_i == 0:
				print("\nAnalysing full spectra")
			else:
				print("\nAnalysing spectra with copy number", s_i-1)
			s.analyse(min_perc=min_perc, min_elem=min_elem, verbose=verbose)
			if s.peaks:
				self.limy = int(max(int(s.maxval * 1.1 / 1000) * 1000, self.limy))
				self.limx = int(max(min(s.peaks[-1].mean * 2, len(s.histogram)), self.limx))

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
		goal = 0.99 * sum([x.elements for x in general_dists])
		maxpeaks = 10
		general_dists.sort(key=lambda x: -x.elements)
		af = []
		peaks = 0
		covered = 0
		for x in general_dists:
			af.append(x.mean)
			peaks += 1
			covered += x.elements
			if peaks == maxpeaks or covered > goal:
				break

		# step 3, report for each peak
		# get the candidate peak on each spectra
		for f in af:
			total = 0
			pd = {}
			for i in range(len(self.spectras) - 1):
				m = [(x.mean, x.elements) for x in self.spectras[1 + i].peaks if x.mean > 0.8 * f and x.mean < 1.2 * f]
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
			opt_freq = int(self.spectras[0].peaks[hpi-1].mean)
			absent_count = self.spectras[1].histogram[opt_freq]
			present_count = self.spectras[2].histogram[opt_freq]
			return (present_count / (absent_count + present_count)) * 100.0
		else:
			return 0.0

def get_properties_from_file(input_file):
	k = 27
	mx = False

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
		i+=1
	f.close()

	return k, mx

def main():
	# ----- command line parsing -----
	parser = argparse.ArgumentParser(
		description="""Analyse a comp matrix file with respect to the distributions and copy numbers seen within.""")

	parser.add_argument("input", type=str,
						help="The input should be either a KAT spectra-cn matrix file or a KAT histogram file.")

	parser.add_argument("-c", "--cns", type=int, default=4,
						help="The number of copy numbers to consider in the analysis.  Only applicable if input is a matrix file.")
	parser.add_argument("-f", "--freq_cutoff", type=int, default=500,
						help="The maximum frequency cutoff point to consider.  Analysis will be done up to this frequency.")
	parser.add_argument("-p", "--min_perc", type=int, default=1,
						help="Any new distribution that adds less to min perc kmers on the iterative analysis will not be added.")
	parser.add_argument("-e", "--min_elem", type=int, default=100000,
						help="Any new distribution that adds less to min elem kmers on the iterative analysis will not be added.")
	parser.add_argument("--plot", action='store_true',
						help="Plot fitted distributions to each peak.")
	parser.add_argument("-z", "--homozygous_peak", type=int, default=0,
						help="The approximate kmer frequency for the homozygous peak.  Allows us to calculate a more accurate genome size estimate.")
	parser.add_argument("-v", "--verbose", action='store_true',
						help="Print additional information.")

	args = parser.parse_args()

	k, mx = get_properties_from_file(args.input)
	print("Input file generated using K", k)
	print("Input is a", "matrix" if mx else "histogram", "file")

	a = None
	if mx:
		print("Processing", args.cns, "spectra")
		a = MXKmerSpectraAnalysis(args.input, cns_cutoff=args.cns, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)
	else:
		a = HistKmerSpectraAnalysis(args.input, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)

	if not a:
		raise ValueError("Couldn't generate a valid spectra analysis")

	a.analyse(min_perc=args.min_perc, min_elem=args.min_elem, verbose=args.verbose)
	a.peak_stats()
	if args.plot:
		a.plot()

if __name__ == '__main__':
	main()