#!/usr/bin/env python3

import argparse
import abc
import sys
import traceback
import time

import numpy as np
from scipy import mean, optimize
from scipy.stats import norm
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


R2PI = np.sqrt(2.0 * np.pi)
PEAK_TOLERANCE = 10.0	# Higher number requires actual peak to be closer to expected peak (i.e. stricter)




class KmerSpectra(object):
	"""
	A kmer spectra, comprised of different peaks. Contains the general fitting method.
	"""
	def __init__(self, histogram, k=27):
		"""
		Inititalise the spectra with the actual histogram to model
		:param histogram: Histogram derived from one of the KAT tools
		:param k: K value used to construct the histogram
		"""

		self.histogram = histogram
		self.realTy = np.array(histogram)
		self.k = k		# K-mer value used.  Comes in handy for genome size estimations and heterzygous rates
		self.peaks = []  # Will contain all distributions that model the histogram
		self.Tx = np.linspace(0, len(self.histogram) - 1, len(self.histogram))
		self.Ty = np.zeros_like(self.Tx)
		self.fmax = 0	# Position of global maxima in actual histogram
		self.fmin = 0	# Position of global minima in actual histogram
		self.maxval = 0	# Value at global maxima (shortcut for self.histogram[self.fmax])

	def analyse(self,min_perc=1, min_elem=100000, verbose=False):
		"""
		Analyse the histogram for peaks
		:param min_perc:
		:param min_elem:
		:param verbose:
		:return:
		"""

		if verbose:	print("Creating initial peaks...", end="", flush=True)
		self.create_peaks()

		if self.peaks:
			if verbose:
				print("done.", len(self.peaks), "peaks initially created:")
				self.printPeaks()
				print("Global minima @ Frequency=" + str(self.fmin) + " (" + str(self.histogram[self.fmin]) + ")")
				print("Global maxima @ Frequency=" + str(self.fmax) + " (" + str(self.histogram[self.fmax]) + ")")
				print("Locally optimising each peak...", end="")
			self.locally_optimise_peaks()
			if verbose:
				print("done.")
				self.printPeaks()
				print("Optimising cumulative distribution to histogram...", end="", flush=True)
			try:
				self.optimise()
				if verbose:
					print("done.")
					self.printPeaks()
			except Exception as inst:
				print("WARNING: problem optimising peaks. It is likely that the spectra is too complex to analyse properly.  Output for this spectra may not be valid.", file=sys.stderr)
				print(inst, file=sys.stderr)
				pass
		elif verbose:
			print("done. No peaks created")


	def smooth(self, x, window_len=3):
		"""
		Smooths the histogram using a moving average
		:param x: Histogram to smooth
		:param window_len: Window length, larger value is smoother.  min (and default) is 3
		:return: A smoothed version of x
		"""
		if x.ndim != 1:
			raise ValueError("Smooth only accepts 1 dimension arrays.")

		if x.size < window_len or window_len < 3:
			return x

		s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
		w = np.ones(window_len, 'd')
		y = np.convolve(w / w.sum(), s, mode='valid')
		return y


	def create_peaks(self):
		"""
		Creates a set of peaks based on the global maxima (ignoring the likely first peak at low k-mer frequency).
		Because we expect this to be a poisson distribution we assume peaks can be found at multiples of the global
		maxima.
		:return:
		"""

		# Walk till first local minimum (d(f)>0)
		# Also double check the following two steps, rather than just the next one.
		# Sometimes we can get a strange laddering affect in alternate frequencies which prevent
		# us from correctly detecting the minima
		fmin = 0
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

		if fmax < 10:
			# Nothing to do in this case.  Not enough data to create peaks.
			return

		# Set global histogram properties to class variables for later use
		self.fmax = fmax
		self.fmin = fmin
		self.maxval = self.histogram[fmax]

		# Unless otherwise specified we assume fmax represents the homozygous peak (primary content)
		# Explore expected peak sites, also look for heterozygous content.
		for mu in [fmax / 3.0, fmax / 2.0, fmax, fmax * 2, fmax * 3, fmax * 4, fmax * 5]:

			# This seems like a reasonable guess for the initial stddev of the peak, given the
			# behaviour of the poisson distribution
			sigma = mu / 6.0

			# We are interested in a region up to 2 stddevs from the mean (95% coverage)
			radius = sigma * 2.0

			# Conditions:
			# - we need at least a radius of 2
			# - f must be greater than fmin
			# - The extent of the distribution including the radius should extend over the histogram limits
			if int(radius) >= 2 and int(mu) > fmin and int(mu - radius) > 0 and int(mu + radius) < len(self.histogram):

				# This code assumes a maxima exists here
				self.peaks.append(KmerPeak(
					int(mu),  			# Mean
					sigma,  			# Sigma
					sum([self.histogram[j] for j in range(int(mu - radius), int(mu + radius))]),  # Actual elements from histogram (not distribution)
					self.histogram[int(mu)],	# Peak value at mean
					mu == fmax			# Whether or not this is the primary peak
				)) 	# Right


				# This code will only return a peak if it can find a maxima
				#lm = self.find_maxima(int(f), int(radius), min_perc, min_elem)
				#if lm:
				#	self.add_peak_and_update_cuts(lm, reset_opt=True)


	# ----------------- OPTIMIZATION --------------------

	def locally_optimise_peaks(self):

		# Make a copy of the histogram and modify...
		base = np.array(self.histogram)

		# Locally optimize each peak,
		# the cumulative fitted distribution of all peaks will likely overshoot the actual histogram
		# at this stage
		for p_i, p in enumerate(self.peaks):
			p.optimise(base)


		# For debugging
		plt.plot(self.histogram, color='black')
		for p_i, p in enumerate(self.peaks):
			plt.plot(p.Ty)
		plt.xlim(0,70)
		plt.ylim(0,120000000)
		plt.show()


	def update_fitted_histogram(self):
		"""
		This function updates the fitted histogram based on the current parameters in each of the peaks in this spectra
		:return: The newly fitted histogram (self.fitted_histogram)
		"""
		# TODO there's probably a super fast numpy vectorised way of doing this.
		self.Ty = np.zeros_like(self.Tx)
		for p in self.peaks:
			scaling_factor = p.peak / gaussian(p.mean, p.mean, p.sigma)
			for i, x in enumerate(self.Tx):
				self.Ty[i] += int(gaussian(x, p.mean, p.sigma) * scaling_factor)

		return self.Ty


	def objective_func(self, p):
		"""
		Our objective is to create a set of distributions that fits the real histogram as closely as possible
		We do this by trying to minimise the difference between our fitted histogram (cumulative sum of
		all distributions) and the real histogram.  The smaller the difference the better.
		:param p: New set of parameters adjusted by the optimiser
		:return: Scalar values representing the
		"""
		if len(p) != len(self.peaks) * 2:
			raise ValueError("Parameters and peaks have got out of sync")

		for i in range(len(self.peaks)):
			self.peaks[i].sigma = p[i * 2]
			self.peaks[i].peak = p[i * 2 + 1]

		# Quick check through each peak to accumulate any penalties
		penalty = 0
		for p in self.peaks:
			if p.left() < 1:
				penalty += (abs(1 - p.left()) * 10000000000.0)
			if p.peak > self.histogram[p.mean]:
				penalty += (abs(p.peak - self.histogram[p.mean]) * 1000000000.0)

		# If we've added a penalty then just return that
		if penalty:
			return penalty

		# Recalculate the fitted histogram based on information in the new parameters "p"
		self.update_fitted_histogram()

		# Create a list of differences between actual and fitted histogram in the area of interest
		delta = np.abs(self.realTy - self.Ty)

		score = sum(delta)

		return score

	def optimise(self):
		p = []
		bounds = []
		for pk in self.peaks:
			p.append(pk.sigma)
			p.append(float(pk.peak))

		# Reset Tx and Ty in case the histogram has been modified
		self.Tx = np.linspace(0, len(self.histogram) - 1, len(self.histogram))
		self.Ty = np.zeros_like(self.Tx)
		self.realTy = np.array(self.histogram)

		# Optimise
		res = optimize.fmin_slsqp(self.objective_func, p, full_output=True)

		if res[-1] < 1 or res[-1] > 4:
			raise RuntimeError("It is likely that the spectra is too complex to analyse properly.  Stopping analysis.\nOptimisation results:\n" + str(res[-2]))

		# once the better fit is found, check if by taking the unfitted elements new distributions arise.
		return


	def getHomozygousPeakIndex(self, approx_freq=0):
		# Work out which peak is the homozygous peak
		peak_index = 3
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

		hom_peak_index = 3 # len(self.peaks) if hom_peak == 0 else hom_peak

		if hom_peak_index == 0:
			return 0

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

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

class KmerPeak(object):
	"""A distribution representing kmers covered a certain number of times.
	Contains methods for fitting to an interval"""

	def __init__(self, mean, sigma, elements, peak, primary):
		self.mean = mean
		self.sigma = sigma
		self.elements = elements
		self.peak = peak
		self.primary = primary

	def left(self):
		return int(self.mean - self.radius())

	def right(self):
		return int(self.mean + self.radius())

	def radius(self):
		"""
		Returns the radius we are interested in for this peak.
		Current set at 2 * the stddev, so 2 * radius gives 95% coverage
		:return: region of interest from mean
		"""
		return 2.0 * self.sigma

	def __str__(self):
		return "Peak of " + str(self.peak) + " at frequency " + str(int(self.mean)) + ", with volume of " + \
			str(int(self.elements)) + " elements between frequencies of " + str(self.left()) + " and " + str(self.right()) + ".  StdDev: " + str(self.sigma) + "; Primary: " + str(self.primary)

	def toTabString(self):
		return "\t".join([str(self.left()), str(int(self.mean)), str(self.right()), str(int(self.peak)), str(int(self.elements)), str(self.primary)])

	@staticmethod
	def getTabHeader():
		return "Left\tMean\tRight\tMax\tVolume\tPrimary"

	def update_fitted_distribution(self):
		"""
		Updates both the distribution and
		:param scaling_factor:
		:return: updated fitted histogram
		"""

		# TODO there's probably a super fast numpy vectorised way of doing this.
		for i, x in enumerate(self.Tx):
			self.Ty[i] = gaussian(x, self.mean, self.sigma) * self.scaling_factor

		self.elements = sum(self.Ty)
		self.scaling_factor = float(self.peak) / gaussian(self.mean, self.mean, self.sigma)

		return self.Ty

	def objective_func(self, p):
		"""
		Fit this gaussian distribution as closely as possible to the histogram using the given parameters
		:param p: The parameters to use
		:return: Value representing the delta between the fitted gaussian and the histogram
		"""
		self.sigma = p[0]

		# If updated sigma causes left() is < 1, then return a huge penalty
		if self.left() < 1:
			return abs(1 - self.left()) * 10000000000.0

		# Updates the histogram represented by this specific peak based on the new sigma
		self.update_fitted_distribution()

		# Return the distance between the fitted peak and the actual histogram at each site
		delta = np.abs(self.Ty - self.histogram)

		# Sum the distances to provide overall level of difference
		return sum(delta)

	def optimise(self, histogram):
		"""
		Tries to fit this single guassian distribution to this point in the histogram as closely as possible
		:param histogram:
		:return:
		"""

		# Sanity check... this should never happen though.
		if len(histogram) == 0:
			raise RuntimeError("Encountered peak with no range")

		# Save histogram for easy access
		self.histogram = histogram

		# Initialise element count to within radius (2 stddev) of the mean from the real histogram
		self.elements = sum([x for x in self.histogram[self.left() : self.right()]])

		# Calculate scaling factor based on std-dev and number of elements in this guassian
		self.scaling_factor = float(self.peak) / gaussian(self.mean, self.mean, self.sigma)

		# Create Tx and Ty (represents fitted histograms)
		self.Tx = np.linspace(0, len(histogram) - 1, len(histogram))
		self.Ty = np.zeros_like(self.Tx)

		# Update the fitted histograms based on this peak's scaled gaussian
		self.update_fitted_distribution()

		# Set up variables to optimise (just sigma)
		p = [float(self.sigma)]

		# Set the optimal sigma and peak values for this
		res = optimize.leastsq(self.objective_func, p, full_output=True, maxfev=100)
		if res[-1] < 1 or res[-1] > 4:
			print("It is likely that the spectra is too complex to analyse properly.  Stopping analysis.\nOptimisation results:\n" + str(res[-2]))

		return

def plot_hist(h, xmax, ymax, label="", xlab='Kmer Frequency', ylab='# Distinct Kmers'):
	plt.plot([min(ymax, x) for x in h[:xmax]], label=label)
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
			plot_hist(self.spectra.histogram, xmax, ymax, label="Histogram")
			plot_hist(self.spectra.update_fitted_histogram(), xmax, ymax, label="Fitted distribution")

			for p_i, p in enumerate(self.spectra.peaks, start=1):
				plot_hist(p.Ty, xmax, ymax, label="fit dist %d" % p_i)

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

			#self.cov_spectra.printPeaks()

			fig = plt.figure()
			plot_hist(self.cov_spectra.histogram, xmax, ymax, label="Histogram")
			plot_hist(self.cov_spectra.update_fitted_histogram(1, xmax + 1), xmax, ymax, label="Fitted distribution")

			for p_i, p in enumerate(self.cov_spectra.peaks, start=1):
				plot_hist(p.points(1, xmax + 1), xmax, ymax, label="fit dist %d" % p_i)

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
			plot_hist(self.gc_dist.histogram, xmax, ymax, label="Histogram", xlab="GC count")
			plot_hist(self.gc_dist.update_fitted_histogram(1, xmax + 1), xmax, ymax, label="Fitted distribution", xlab="GC count")

			for p_i, p in enumerate(self.gc_dist.peaks, start=1):
				plot_hist(p.points(1, xmax + 1), xmax, ymax, label="fit dist %d" % p_i, xlab="GC count")

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

			plot_hist(s.histogram, xmax, ymax, label=slabel)
			plot_hist(s.update_fitted_histogram(1, xmax + 1), xmax, ymax, label=slabel + " fit")
			for p_i, p in enumerate(s.peaks, start=1):
				plot_hist(p.points(1, xmax + 1), xmax, ymax, label="fit dist %d" % p_i)

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
