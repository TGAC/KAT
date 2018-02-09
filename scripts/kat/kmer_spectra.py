import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

from .kmer_peak import KmerPeak


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
		self.k = k  # K-mer value used.  Comes in handy for genome size estimations and heterzygous rates
		self.peaks = []  # Will contain all distributions that model the histogram
		self.Tx = np.linspace(0, len(self.histogram) - 1, len(self.histogram))
		self.Ty = np.zeros_like(self.Tx)
		self.fmax = 0  # Position of global maxima in actual histogram
		self.fmin = 0  # Position of global minima in actual histogram
		self.maxval = 0  # Value at global maxima (shortcut for self.histogram[self.fmax])

	def analyse(self, min_perc=1, min_elem=100000, verbose=False):
		"""
		Analyse the histogram for peaks
		:param min_perc:
		:param min_elem:
		:param verbose:
		:return:
		"""

		if verbose:
			print()
			print("Creating initial peaks...", end="", flush=True)
		self.create_peaks()

		if self.peaks:
			if verbose:
				print("done.", len(self.peaks), "peaks initially created")
				print()
				self.printPeaks()
				print()
				print("Global minima @ Frequency=" + str(self.fmin) + " (" + str(self.histogram[self.fmin]) + ")")
				print("Global maxima @ Frequency=" + str(self.fmax) + " (" + str(self.histogram[self.fmax]) + ")")
				print()
				print("Locally optimising each peak...", end="")
			for p_i, p in enumerate(self.peaks):
				p.optimise(self.histogram)
				# For debugging
			if True:
				plt.plot(self.histogram, color='black')
				for p_i, p in enumerate(self.peaks):
					plt.plot(p.Ty)
				plt.xlim(0, 70)
				plt.ylim(0, 120000000)
				plt.show()

			if verbose:
				print("done.")
				print()
				self.printPeaks()
				print()
				print("Optimising cumulative distribution to histogram...", end="", flush=True)
			try:
				self.optimise()
				if verbose:
					print("done.")
					print()
					self.printPeaks()
			except Exception as inst:
				print(
					"WARNING: problem optimising peaks. It is likely that the spectra is too complex to analyse properly.  Output for this spectra may not be valid.",
					file=sys.stderr)
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
			if self.histogram[i] < self.histogram[i + 1] and self.histogram[i] < self.histogram[i + 2]:
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
		for mu in [fmax / 2.0, fmax, fmax * 2, fmax * 3, fmax * 4, fmax * 5]:

			# In a poisson distribution the mean is the variance, so the stddev is simply the square root of the mean
			sigma = np.sqrt(mu)

			# We are (at present, only) interested in a region up to 2 stddevs from the mean (95% coverage)
			radius = int(sigma * 2.0)
			mean = int(mu)

			# Conditions:
			# - we need at least a radius of 2
			# - f must be greater than fmin
			# - The extent of the distribution including the radius should extend over the histogram limits
			if radius >= 2 and mean > fmin and mu - radius > 0 and mu + radius < len(self.histogram):
				# This code assumes a maxima exists here
				self.peaks.append(KmerPeak(
					mean,		# Mean
					sigma,  	# Std dev
					self.histogram[mean], # Use the histogram value at this position as an initial value for this peak
					mean == fmax  # Whether or not this is the primary peak
				))

			# This code will only return a peak if it can find a maxima
			# lm = self.find_maxima(int(f), int(radius), min_perc, min_elem)
			# if lm:
			#	self.add_peak_and_update_cuts(lm, reset_opt=True)

	def _updateModel(self, params):
		"""
		This function updates the fitted histogram based on the current parameters in each of the peaks in this spectra
		:return: The newly fitted histogram (self.fitted_histogram)
		"""
		# TODO there's probably a super fast numpy vectorised way of doing this.
		self.Ty = np.zeros_like(self.Tx)
		for i in range(len(self.peaks)):
			self.Ty += self.peaks[i].updateModel(params[i * 2], params[i * 2 + 1])

		return

	def _objectiveFunc(self, params):
		"""
		Our objective is to create a set of distributions that fits the real histogram as closely as possible
		We do this by trying to minimise the difference between our fitted histogram (cumulative sum of
		all distributions) and the real histogram.  The smaller the difference the better.
		:param params: New set of parameters adjusted by the optimiser
		:return: Scalar values representing the
		"""

		# Quick sanity check (probably can drop this to save time)
		if len(params) != len(self.peaks) * 2:
			raise ValueError("Parameters and peaks have got out of sync")

		# Recalculate the fitted histogram based on information in the new parameters
		self._updateModel(params)

		# Create a list of differences between actual and fitted histogram in the area of interest
		delta = self.Ty - self.realTy

		# We want to heavily penalise all points which exceed the histogram more harshly than those underneath
		for i in range(len(delta)):
			d = delta[i]
			if d > 0:
				delta[i] = d * 2#np.power(d, 2)
			elif d < 0:
				delta[i] = d

		score = sum(delta)

		return score

	def optimise(self):
		"""
		Given the full set of peaks, adjust all their heights in order to best fit the acutal histogram
		We also put some bounds around the limits these values can take to stop them going crazy
		"""
		params = []
		bounds = []
		for p in self.peaks:
			params.append(p.peak())
			bounds.append((0.1 * self.histogram[p.mean()], self.histogram[p.mean()]))
			params.append(p.stddev())
			bounds.append((np.sqrt(p.mean()), p.mean()))


		# Reset Tx and Ty in case the histogram has been modified
		self.Tx = np.linspace(0, len(self.histogram) - 1, len(self.histogram))
		self.Ty = np.zeros_like(self.Tx)
		self.realTy = np.array(self.histogram)

		# Optimise
		res = optimize.minimize(self._objectiveFunc, np.array(params), bounds=bounds, method="L-BFGS-B")
		if not res.success:
			raise RuntimeError(
				"It is likely that the spectra is too complex to analyse properly.  Stopping analysis.\nOptimisation results:\n" + str(
					res[-2]))

		# once the better fit is found, check if by taking the unfitted elements new distributions arise.
		return

	def getHomozygousPeakIndex(self, approx_freq=0):
		"""
		If an approximate frequency is not provided then we assume the largest peak is the homozygous peak
		:param approx_freq: User provided guide for roughly where the homozygous peak should be located
		:return: The 1-based index of the primary peak
		"""
		# Work out which peak is the homozygous peak
		if approx_freq > 0:
			min_peak_index = 0
			delta_peak_freq = 1000000
			for p_i, p in enumerate(self.peaks, start=1):
				delta = abs(p.mean() - approx_freq)
				if delta_peak_freq > delta:
					delta_peak_freq = delta
					min_peak_index = p_i
			return min_peak_index
		else:
			for i, p in enumerate(self.peaks, start=1):
				if p.primary:
					return i

		return 0

	def calcGenomeSize(self, hom_peak=0):
		"""
		Attempts to calculate the genome size.  Requires knowledge of which peak represents the homozygous peak to
		work.  Essentially we sum the volume under the heterzygous and homozygous peaks, then multiply the volume under
		peaks representing repeat content by the relative index after the homozygous peak.
		:param hom_peak: User provided guide for roughly where the homozygous peak should be located, if 0, then we assume it's the largest peak
		:return: The estimated genome size
		"""

		hom_peak_index = self.getHomozygousPeakIndex(hom_peak)

		if hom_peak_index == 0:
			return 0

		sum = 0
		for p_i, p in enumerate(self.peaks, start=1):
			if p_i > hom_peak_index:
				sum += (p_i - hom_peak_index) * p.elements()
			else:
				sum += p.elements()

		return sum

	def calcHetRate(self, genome_size=0, hom_peak=0):
		"""
		Calculate the heterozygous rate based on the fraction of the whole genome falling into the heterozygous peak
		:param genome_size: User provided genome size, if 0 we try to calculate it ourselves
		:param hom_peak: User provided guide for roughly where the homozygous peak should be located, if 0, then we assume it's the largest peak
		:return: The heterozygous rate
		"""

		genomesize = genome_size if genome_size > 0 else self.calcGenomeSize()
		hom_peak_index = self.getHomozygousPeakIndex(hom_peak)

		# First do a sanity check to make sure there is some heterzygous content to work with
		if hom_peak_index < 2:
			return 0.0

		sum = 0
		for p_i, p in enumerate(self.peaks, start=1):
			# Skip the last peak
			if p_i >= hom_peak_index:
				break
			sum += p.elements() / self.k

		return (sum / genomesize) * 100.0

	def calcKmerCoverage(self):
		tot_vol = sum([x.elements() for x in self.peaks])
		weighted = sum([x.mean() * x.elements() for x in self.peaks])
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
