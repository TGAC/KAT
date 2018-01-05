#!/usr/bin/env python3

import argparse
import abc
import sys
import traceback
import time

import numpy as np
from scipy import mean, optimize
import matplotlib.pyplot as plt


R2PI = np.sqrt(2.0 * np.pi)

class KmerSpectra(object):
	"""A kmer spectra, comprised of different peaks.
	Contains the general fitting method"""

	###----------------- INIT, LOAD, DUMP, ETC --------------------
	def __init__(self, histogram, k=27):
		"""Init with the objective histogram"""
		self.histogram = histogram
		self.peaks = []
		self.k = k
		self.fitted_histogram = [0] * len(self.histogram)
		self.fmax = 0
		self.fmin = 0
		self.maxval = 0

	def total_values(self, start=1, end=10000):
		self.fitted_histogram = list(map(sum, list(zip(*[x.points(start, min(end, len(self.histogram))) for x in self.peaks]))))
		return self.fitted_histogram

	###----------------- FIND AND CREATE PEAKS --------------------
	def analyse(self,min_perc=1, min_elem=100000, verbose=False, gcd=False):
		if verbose:
			print("Creating initial peaks...", end="")
			sys.stdout.flush()
		self.create_peaks(min_perc=min_perc, min_elem=min_elem, verbose=verbose, gcd=gcd)

		if self.peaks:
			if verbose:
				print("done.", len(self.peaks), "peaks initially created:")
				self.printPeaks()
				print("Global minima @ Frequency=" + str(self.fmin) + " (" + str(self.histogram[self.fmin]) + ")")
				print("Global maxima @ Frequency=" + str(self.fmax) + " (" + str(self.histogram[self.fmax]) + ")")
				print("Optimising peaks...", end="")
			self.optimize_peaks(min_perc=min_perc, min_elem=min_elem, verbose=verbose)
			if verbose:
				print("done. There are now", len(self.peaks), "peaks:")
				self.printPeaks()
				print("Optimising all peaks in spectra...", end="")
				sys.stdout.flush()
			try:
				self.optimize_overall()
				if verbose:
					print("done.", len(self.peaks), "peaks present after optimisation:")
					self.printPeaks()
			except:
				print("WARNING: problem optimising peaks. It is likely that the spectra is either too complex to analyse properly.  Output for this spectra may not be valid.", file=sys.stderr)
				pass
		elif verbose:
			print("done. No peaks created")

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
			print("Distribution at %d too small to be considered (%d elements)" % (start + fmax, sum(hs)))
			return 0
		# TODO: do further validation
		# print "maxima found on %d, for %s" % (fmax,hs)
		return start + fmax

	def add_peak_and_update_cuts(self, lm, reset_opt=False):
		# Receives a local maxima, adds a peak there if there was none, updates the cuts.
		# If reset_opt is set, it will reset all the peak counts/sd and only keep the means

		# Check if peak is below fmin.  If so skip it.
		if lm < self.fmin:
			return False

		# Validation.  Extract means for all peaks and ensure they are not too close to one another
		fdists = [x.mean for x in self.peaks]
		for f in fdists:
			if lm >= f - f / 4 and lm <= f + f / 4:
				#print("WARNING!!! Dist at %d is not to be added, because existing dist at %d is too close to it." % (lm, f))
				return False

		# Add the new peak and sort
		fdists.append(lm)
		fdists.sort()

		# Use what is hopefully the cutoff frequency as the start of the first peak
		# Then try to segment the spectra into different sections based on the peak frequencies
		self.cuts = [self.fmin]
		for i in range(len(fdists) - 1):
			self.cuts.append(int(fdists[i] * (1 + float(fdists[i]) / (fdists[i] + fdists[i + 1]))))
		self.cuts.append(int(min(len(self.histogram) - 1, fdists[-1] * 1.5)))
		self.cuts.sort()
		#print("Cuts: %s" % str(self.cuts))
		if reset_opt:
			self.peaks = []

		for i in range(len(fdists)):
			if reset_opt or i == fdists.index(lm):
					self.peaks.append(KmerPeak(
											fdists[i], 			# Mean
											fdists[i] / 6,		# Shape
											sum([self.histogram[j] for j in range(self.cuts[i], self.cuts[i+1])]),	# Elements
											self.histogram[fdists[i]], #Peak
											self.cuts[i], 		# Left
											self.cuts[i+1]))	# Right
					self.peaks.sort(key=lambda x: x.mean)
		return True

	def create_peaks(self, min_perc, min_elem, verbose=False, gcd=False):

		# walk till first local minimum (d(f)>0)
		# Double check the following two steps, rather than just the next one.
		# Sometimes we can get a strange laddering affect in alternate frequencies which prevent
		# us from correctly detecting the minima
		fmin = 0
		if not gcd:
			for i in range(1, len(self.histogram) - 2):
				if self.histogram[i] < self.histogram[i+1] and self.histogram[i] < self.histogram[i+2]:
					fmin = i
					break

		# Sometimes we might not find a local minima, it depends on what sort of content is in the spectra.
		# In this case just reset all spectra measures
		if not fmin:
			self.fmax = 0
			self.maxval = 0

		# Find the global fm|p(fm)=max(p)
		fmax = self.histogram.index(max(self.histogram[fmin:]))

		if not gcd and fmax < 10:
			#print("Kmer count at max frequency  is < 10.  Not enough data to perform analysis.")
			return

		# Set measures to class variables
		self.fmax = fmax
		self.fmin = fmin
		self.maxval = self.histogram[fmax]

		# Explore fm/x and fm*x for x in [1,4]
		for f in [fmax / 4.0, fmax / 3.0, fmax / 2.0, fmax, fmax * 2, fmax * 3, fmax * 4]:
			radius = f / 3.0
			if int(radius) > 0 and int(f - radius) > 0 and int(f + radius) < len(self.histogram):
				lm = self.find_maxima(int(f), int(radius), min_perc, min_elem)
				if lm:
					self.add_peak_and_update_cuts(lm, reset_opt=True)

	# ----------------- OPTIMIZATION --------------------

	def optimize_peaks(self, min_perc, min_elem, verbose=False):
		sortedpeaks = [x for x in self.peaks]
		sortedpeaks.sort(key=lambda x: -x.elements)

		# create a base as the histogram and start from there
		base = [x for x in self.histogram]
		for p_i, p in enumerate(sortedpeaks):
			# locally optimize the peak
			p.constrained_opt(p.left, base[p.left:p.right])
			# substract the peak effect from the baseline
			for i in range(len(self.histogram)):
				base[i] = base[i] - p.point(i)

		updated = False
		for f in [self.fmax / 4, self.fmax / 3, self.fmax / 2, self.fmax, self.fmax * 2, self.fmax * 3, self.fmax * 4]:
			if f <= 5:
				continue

			if int(f + f / 5) < len(self.histogram) and sum(
					[base[x] for x in range(int(f - f / 5), int(f + f / 5))]) > .01 * sum(
					[p.elements for p in self.peaks]):
				lm = self.find_maxima(int(f), int(f / 5), min_perc, min_elem, base)
				if lm:
					if self.add_peak_and_update_cuts(lm, reset_opt=False):
						updated = True
		if updated: self.optimize_peaks(min_perc, min_elem)

	def residues(self, p):
		"""p has a triplet for each distribution"""
		# print p
		for i in range(len(self.peaks)):
			self.peaks[i].mean = p[i * 5]
			self.peaks[i].shape = p[i * 5 + 1]
			self.peaks[i].elements = p[i * 5 + 2]
			self.peaks[i].left = p[i * 5 + 3]
			self.peaks[i].right = p[i * 5 + 4]

		self.total_values()
		r = [self.fitted_histogram[i] - self.histogram[i] for i in range(self.cuts[0], self.cuts[-1])]
		return r

	def optimize_overall(self):
		p = []
		for pk in self.peaks:
			p.append(pk.mean)
			p.append(pk.shape)
			p.append(pk.elements)
			p.append(pk.left)
			p.append(pk.right)

		# Optimise
		res = optimize.leastsq(self.residues, p, full_output=True)
		if res[-1] < 1 or res[-1] > 4:
			raise RuntimeError("It is likely that the spectra is too complex to analyse properly.  Stopping analysis.\nOptimisation results:\n" + str(res[-2]))

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

		if hom_peak_index == 0:
			return 0

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
		return int(weighted / tot_vol) if tot_vol > 0 else 0

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
		print("Homozygous peak index:", hp, ("" if hom_peak_freq > 0 else "(assumed)"))
		gs = self.calcGenomeSize(hom_peak=hp)
		print("Estimated genome size:", '{0:.2f}'.format(float(gs) / 1000000.0), "Mbp")
		if (hp > 1):
			hr = self.calcHetRate(gs)
			print("Estimated heterozygous rate:", "{0:.2f}".format(hr) + "%")


class KmerPeak(object):
	"""A distribution representing kmers covered a certain number of times.
	Contains methods for fitting to an interval"""

	def __init__(self, mean, shape, elements, peak, left, right):
		self.mean = mean
		self.shape = shape
		self.elements = elements
		self.peak = peak
		self.left = left
		self.right = right

	def __str__(self):
		return "Peak of " + str(self.peak) + " at frequency " + str(int(self.mean)) + ", with volume of " + \
			str(int(self.elements)) + " elements between frequencies of " + str(self.left) + " and " + str(self.right) + ".  Shape: " + str(self.shape)

	def toTabString(self):
		return "\t".join([str(int(self.left)), str(int(self.mean)), str(int(self.right)), str(int(self.peak)), str(int(self.elements))])

	@staticmethod
	def getTabHeader():
		return "Left\tMean\tRight\tMax\tVolume"

	def point(self, x):
		"""Normalized Gaussian"""
		#return float(self.elements) / R2PI / self.shape * np.exp(-(x - self.mean) ** 2 / 2. / self.shape ** 2)
		return self.__p1 * np.exp(-(x - self.mean) ** 2 / 2. / self.shape ** 2)

	def points(self, start, end):
		self.__p1 = float(self.elements) / R2PI / self.shape
		return [self.point(x) for x in range(start, end)]

	def residues(self, p, offset, histogram):
		"""residues of fitting self to offset obj[0] with values obj[1]"""
		self.__p1 = float(self.elements) / R2PI / self.shape
		r = [self.point(offset + i) - histogram[i] for i in range(len(histogram))]
		return r

	def constrained_opt(self, offset, histogram):
		"""Does a constrained optimization of fitting to match the histogram"""
		if len(histogram) == 0:
			raise RuntimeError("Encountered peak with no range")
		optimize.leastsq(self.residues, (self.mean, self.shape, self.elements), (offset, histogram))
		pass

def plot_hist(h, points, cap, label="", to_screen=False, to_file=None, xlab='Kmer Frequency', ylab='# Distinct Kmers'):
	plt.plot([min(cap, x) for x in h[:points]], label=label)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.legend()



def plot_hist_df(self, h, points, cap, xlab='Kmer Frequency', ylab='# Distinct Kmers'):
	plt.plot([max(-cap, min(cap, x)) for x in [h[i + 1] - h[i] for i in range(points)]])
	plt.xlabel(xlab)
	plt.ylabel(ylab)
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
	def plot(self, points=0, cap=0, to_screen=False, to_files=None):
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

	def plot(self, points=0, cap=0, to_screen=False, to_files=None):
		if 0 == points: points = self.limx
		if 0 == cap: cap = self.limy
		print()
		print("Creating plots")
		print("--------------")
		print()

		if len(self.cov_spectra.peaks) == 0:
			print("No peaks in K-mer frequency histogram.  Not plotting.")
		else:
			print("Plotting K-mer frequency distributions...")
			#self.spectra.printPeaks()

			fig = plt.figure()
			plot_hist(self.spectra.histogram, points, cap, label="Histogram")
			plot_hist(self.spectra.total_values(1, points + 1), points, cap, label="Fitted distribution")

			for p_i, p in enumerate(self.spectra.peaks, start=1):
				plot_hist(p.points(1, points + 1), points, cap, label="fit dist %d" % p_i)

			if to_screen:
				plt.show()

			if to_files:
				filename = to_files + ".dists.png"
				fig.savefig(filename)
				print("- Saved plot to:", filename)

		print()

	def analyse(self, min_perc=1, min_elem=100000, verbose=False):

		# Create the peaks
		if verbose:
			print("Analysing spectra")
		self.spectra.analyse(min_perc=min_perc, min_elem=min_elem, verbose=verbose)
		if self.spectra.peaks:
			self.limy = int(max(int(self.spectra.maxval * 1.1 / 1000) * 1000, self.limy))
			self.limx = int(max(min(self.spectra.peaks[-1].mean * 2, len(self.spectra.histogram)), self.limx))

		if verbose:
			print("Plot limits: y->%d, x->%d" % (self.limy, self.limx))

	def peak_stats(self):
		print()
		print("K-mer frequency spectra statistics")
		print("----------------------------------")
		self.spectra.printGenomeStats(self.hom_peak)

class GCKmerSpectraAnalysis(SpectraAnalysis):
	def __init__(self, filename, freq_cutoff=10000, hom_peak_freq=0, k=27):
		SpectraAnalysis.__init__(self, freq_cutoff=freq_cutoff, hom_peak_freq=hom_peak_freq, k=k)
		cov_histo, gc_histo = self.read_file(filename, freq_cutoff)
		self.mean_gc = sum([i * x for i, x in enumerate(gc_histo)]) / sum(gc_histo)
		self.cov_spectra = KmerSpectra(cov_histo, k=k)
		self.gc_dist = KmerSpectra(gc_histo, k=k)		# Not really a Kmer spectra but let's shoehorn it in anyway


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


	def plot(self, points=0, cap=0, to_screen=False, to_files=None):
		if 0 == points: points = self.limx
		if 0 == cap: cap = self.limy

		print()
		print("Creating plots")
		print("--------------")
		print()
		if len(self.cov_spectra.peaks) == 0:
			print("No peaks in K-mer frequency histogram.  Not plotting.")
		else:

			print("Plotting K-mer frequency distributions...")

			#self.cov_spectra.printPeaks()

			fig = plt.figure()
			plot_hist(self.cov_spectra.histogram, points, cap, label="Histogram")
			plot_hist(self.cov_spectra.total_values(1, points + 1), points, cap, label="Fitted distribution")

			for p_i, p in enumerate(self.cov_spectra.peaks, start=1):
				plot_hist(p.points(1, points + 1), points, cap, label="fit dist %d" % p_i)

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

			points = self.gc_dist.k
			cap = max(self.gc_dist.histogram) * 1.1

			fig = plt.figure()
			plot_hist(self.gc_dist.histogram, points, cap, label="Histogram", xlab="GC count")
			plot_hist(self.gc_dist.total_values(1, points + 1), points, cap, label="Fitted distribution", xlab="GC count")

			for p_i, p in enumerate(self.gc_dist.peaks, start=1):
				plot_hist(p.points(1, points + 1), points, cap, label="fit dist %d" % p_i, xlab="GC count")

			if to_screen:
				plt.show()

			if to_files:
				filename = to_files + ".gc.png"
				fig.savefig(filename)
				print(" - Saved plot to:", filename)

		print()


	def analyse(self, min_perc=1, min_elem=100000, verbose=False):

		# Create the peaks
		if verbose:
			print("Analysing K-mer spectra")
		self.cov_spectra.analyse(min_perc=min_perc, min_elem=min_elem, verbose=verbose)
		if self.cov_spectra.peaks:
			self.limy = int(max(int(self.cov_spectra.maxval * 1.1 / 1000) * 1000, self.limy))
			self.limx = int(max(min(self.cov_spectra.peaks[-1].mean * 2, len(self.cov_spectra.histogram)), self.limx))

		if verbose:
			print("Plot limits: y->%d, x->%d" % (self.limy, self.limx))

		# Create the peaks
		if verbose:
			print("Analysing GC distribution")
		self.gc_dist.analyse(min_perc=1, min_elem=len(self.gc_dist.histogram), verbose=verbose, gcd=True)


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
		self.gc_dist.printPeaks()
		print("Mean GC:", "{0:.2f}".format(self.mean_gc * (100.0 / self.gc_dist.k)) + "%")




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

	def plot(self, points=0, cap=0, to_screen=False, to_files=None):
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
			fig = plt.figure()

			print("Plotting:",slabel)
			s.printPeaks()

			plot_hist(s.histogram, points, cap, label=slabel)
			plot_hist(s.total_values(1, points + 1), points, cap, label=slabel + " fit")
			for p_i, p in enumerate(s.peaks, start=1):
				plot_hist(p.points(1, points + 1), points, cap, label="fit dist %d" % p_i)

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
			elif s_i == 0:
				raise RuntimeError("No peaks detected for full spectra.  Can't continue.")

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
	parser.add_argument("-p", "--min_perc", type=float, default=1.0,
						help="Any new distribution that adds less to min perc kmers on the iterative analysis will not be added.")
	parser.add_argument("-e", "--min_elem", type=int, default=100000,
						help="Any new distribution that adds less to min elem kmers on the iterative analysis will not be added.")
	parser.add_argument("--plot", action='store_true',
						help="Plot fitted distributions to each peak to the screen.")
	parser.add_argument("-z", "--homozygous_peak", type=int, default=0,
						help="The approximate kmer frequency for the homozygous peak.  Allows us to calculate a more accurate genome size estimate.")
	parser.add_argument("-v", "--verbose", action='store_true',
						help="Print additional information.")

	args = parser.parse_args()

	if args.verbose:
		print("\n\nAnalysing distributions for:", args.input)

	k, mx, gcp = get_properties_from_file(args.input)
	if args.verbose:
		print("Input file generated using K", k)

	a = None
	if mx:
		if gcp:
			if args.verbose:
				print("GC vs Coverage matrix file detected")
			a = GCKmerSpectraAnalysis(args.input, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)
		else:
			if args.verbose:
				print("Copy number spectra matrix file detected")
				print("Processing", args.cns, "spectra")
			a = MXKmerSpectraAnalysis(args.input, cns_cutoff=args.cns, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)
	else:
		if args.verbose:
			print("Kmer coverage histogram file detected")
		a = HistKmerSpectraAnalysis(args.input, freq_cutoff=args.freq_cutoff, hom_peak_freq=args.homozygous_peak, k=k)

	if not a:
		raise RuntimeError("Couldn't generate a valid spectra analysis")

	try:
		start = time.time()
		a.analyse(min_perc=args.min_perc, min_elem=args.min_elem, verbose=args.verbose)
		end = time.time()
		print(("" if args.verbose else "done.  ") + "Time taken: ", '{0:.1f}'.format(end - start) + 's')
		a.peak_stats()
		if args.plot or args.output_prefix:
			a.plot(to_screen=args.plot, to_files=args.output_prefix)
	except Exception:
		print("\nERROR\n-----", file=sys.stderr)
		traceback.print_exc(file=sys.stderr)

if __name__ == '__main__':
	main()
