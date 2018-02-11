import abc
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.signal import argrelextrema

from peak import Peak


def smooth(x, window_len=3):
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



class Spectra(object):
	__metaclass__ = abc.ABCMeta

	def __init__(self, histogram, k=27):
		self.histogram = np.array(histogram)
		self.k = k
		self.peaks = None
		self.Tx = np.linspace(0, len(histogram) - 1, len(histogram))
		self.Ty = np.zeros_like(self.Tx)

	@abc.abstractmethod
	def _createInitialPeaks(self, min_perc=1, min_elem=100000, verbose=False):
		pass

	def _updateModel(self, params):
		"""
		This function updates the fitted histogram based on the current parameters in each of the peaks in this spectra
		:return: The newly fitted histogram (self.fitted_histogram)
		"""
		# TODO there's probably a super fast numpy vectorised way of doing this.
		self.Ty = np.zeros_like(self.Tx)
		for i in range(len(self.peaks)):
			if len(params == 2):
				new_mean = self.peaks[i].mean()
				new_peak = params[i * 2]
				new_stddev = params[i * 2 + 1]
			elif len(params == 3):
				new_mean = params[i * 3]
				new_peak = params[i * 3 + 1]
				new_stddev = params[i * 3 + 2]

			self.Ty += self.peaks[i].updateModel(new_mean, new_peak, new_stddev)

		return

	def _residuals(self, params):
		"""
		Our objective is to create a set of distributions that fits the real histogram as closely as possible
		We do this by trying to minimise the difference between our fitted histogram (cumulative sum of
		all distributions) and the real histogram.  The smaller the difference the better.
		:param params: New set of parameters adjusted by the optimiser
		:return: A numpy array of scalar values representing the difference between the model and reality at each X value
		"""

		# Quick sanity check (probably can drop this to save time)
		if len(params) != len(self.peaks) * 2:
			raise ValueError("Parameters and peaks have got out of sync")

		# Recalculate the fitted histogram based on information in the new parameters
		self._updateModel(params)

		# Create a list of differences between actual and fitted histogram in the area of interest
		residuals = self.realTy - self.Ty

		# We want to heavily penalise all points which exceed the histogram more harshly than those underneath
		# But not more harshly than the peak in order to provide a bit more flexiblity in the overall fit.
		# The current multiple was arrived at from trial and error.  Possibly there is a better value to use...
		for i in range(len(residuals)):
			d = residuals[i]
			if d < 0:
				residuals[i] = d * 10

		return residuals

	def optimise(self, vary_mean=False):
		"""
		Given the full set of peaks, adjust all their heights in order to best fit the acutal histogram
		We also put some bounds around the limits these values can take to stop them going crazy
		"""
		if not self.peaks:
			raise ValueError("Can't optimise peaks because none are defined.")

		params = []
		for p in self.peaks:
			if vary_mean:
				params.append(p.mean())
			params.append(p.peak())
			params.append(p.stddev())

		# Reset Tx and Ty in case the histogram has been modified
		self.Tx = np.linspace(0, len(self.histogram) - 1, len(self.histogram))
		self.Ty = np.zeros_like(self.Tx)
		self.realTy = np.array(self.histogram)

		# Optimise
		res = optimize.leastsq(self._residuals, np.array(params), full_output=True)
		if res[-1] < 1 or res[-1] > 4:
			raise RuntimeError(
				"It is likely that the spectra is too complex to analyse properly.  Stopping analysis.\nOptimisation results:\n" + str(
					res[-2]))

		# once the better fit is found, check if by taking the unfitted elements new distributions arise.
		return


	def analyse(self, verbose=False, plot_initial=False):
		"""
		Analyse the histogram for peaks
		:param verbose: Prints additional information about progress to stdout
		"""

		if verbose:
			print()
			print("Creating initial peaks ...", end="", flush=True)
		self._createInitialPeaks()

		if self.peaks:
			if verbose:
				print("done.", len(self.peaks), "peaks initially created")
				print()
				print("Locally optimising each peak ...", end="")
			for p_i, p in enumerate(self.peaks):
				p.optimise(self.histogram)

			# For debugging
			if plot_initial:
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
				print("Fitting cumulative distribution to histogram by adjusting peaks ...", end="", flush=True)
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

	def printPeaks(self):
		if len(self.peaks) > 0:
			print("Index\t" + Peak.getTabHeader())
			for p_i, p in enumerate(self.peaks, start=1):
				print(str(p_i) + "\t" + p.toTabString())
		else:
			print("No peaks detected")





class KmerSpectra(Spectra):
	"""
	A kmer spectra, comprised of different peaks. Contains the general fitting method.
	"""

	def __init__(self, histogram, k=27):
		"""
		Inititalise the spectra with the actual histogram to model
		:param histogram: Histogram derived from one of the KAT tools
		:param k: K value used to construct the histogram
		"""

		# Initialise super
		Spectra.__init__(self, histogram, k)

		# Extra properties for K-mer spectra
		self.fmax = 0  		# Position of global maxima in actual histogram
		self.fmin = 0  		# Position of global minima in actual histogram

	def maxValue(self):
		return self.histogram[self.fmax]



	def _createInitialPeaks(self):
		"""
		Creates a set of peaks based on the global maxima (ignoring the likely first peak at low k-mer frequency).
		Because we expect this to be a poisson distribution we assume peaks can be found at multiples, or fractions, of
		the global maxima.
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
			fmax = 0
		else:
			# Find the global fm|p(fm)=max(p)
			fmax = np.argmax(self.histogram[fmin:,]) + fmin

		# Set member variables
		self.fmin = fmin
		self.fmax = fmax

		if fmax < 10:

			# Not enough data to create peaks in this case.
			self.peaks = None

		else:

			peaks = []

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
					peaks.append(Peak(
						mean,		# Mean
						sigma,  	# Std dev
						self.histogram[mean], # Use the histogram value at this position as an initial value for this peak
						mean == fmax  # Whether or not this is the primary peak
					))

			self.peaks = peaks

		return



	def getHomozygousPeakIndex(self, approx_freq=0):
		"""
		If an approximate frequency is not provided then we assume the largest peak is the homozygous peak
		:param approx_freq: User provided guide for roughly where the homozygous peak should be located
		:return: The 1-based index of the primary peak
		"""
		if approx_freq > 0:
			# User specified a particular frequency to look at, work out which peak is closest and label
			# that the homozygous peak
			min_peak_index = 0
			delta_peak_freq = 1000000
			for p_i, p in enumerate(self.peaks, start=1):
				delta = abs(p.mean() - approx_freq)
				if delta_peak_freq > delta:
					delta_peak_freq = delta
					min_peak_index = p_i
			return min_peak_index
		else:
			# No frequency given.  Use the primary peak (i.e. the one that represents the global maxima)
			for i, p in enumerate(self.peaks, start=1):
				if p.mean() == self.fmax:
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

		hom_peak_index = self.getHomozygousPeakIndex(hom_peak) if hom_peak == 0 else hom_peak

		if hom_peak_index == 0:
			return 0

		sum = 0
		for p_i, p in enumerate(self.peaks, start=1):
			if p_i > hom_peak_index:
				delta = p_i - hom_peak_index + 1
				sum += delta * p.elements()
				p.description = str(delta) + "X"
			elif p_i < hom_peak_index:
				delta = hom_peak_index - p_i + 1
				sum += p.elements() / delta
				p.description = "1/" + str(delta) + "X"
				if p.description == "1/2X":
					p.description = "1/2X (Heterozygous)"
			else:
				sum += p.elements()
				p.description = "1X (Homozygous)"

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


	def printGenomeStats(self, hom_peak_freq):
		hp = self.getHomozygousPeakIndex(hom_peak_freq)

		# This also updates the peaks with information that can be used in labels
		gs = self.calcGenomeSize(hom_peak=hp)

		print("K-value used:", self.k)
		print("Peaks in analysis:", len(self.peaks))
		self.printPeaks()

		print("Mean k-mer frequency:", str(self.calcKmerCoverage()) + "x")
		print("Homozygous peak index:", hp, ("" if hom_peak_freq > 0 else "(assumed)"))
		print("Estimated genome size:", '{0:.2f}'.format(float(gs) / 1000000.0), "Mbp")
		if (hp > 1):
			hr = self.calcHetRate(gs)
			print("Estimated heterozygous rate:", "{0:.2f}".format(hr) + "%")

	def printKmerSpectraProperties(self):
		print("Global minima @ Frequency=" + str(self.fmin) + " (" + str(self.histogram[self.fmin]) + ")")
		print("Global maxima @ Frequency=" + str(self.fmax) + " (" + str(self.histogram[self.fmax]) + ")")


class GCSpectra(Spectra):
	"""
	A kmer spectra, comprised of different peaks. Contains the general fitting method.
	"""

	def __init__(self, histogram, k=27):
		"""
		Inititalise the spectra with the actual histogram to model
		:param histogram: Histogram derived from one of the KAT tools
		:param k: K value used to construct the histogram
		"""

		# Initialise super
		Spectra.__init__(self, histogram, k)


	def _createInitialPeaks(self):
		"""
		Creates a set of peaks based on all maxima after taking the moving average of the histogram
		"""

		#TODO May need to compensate for reduction in elements from moving average
		wlen=3
		smooth_histo = smooth(self.histogram, window_len=wlen)

		# for local maxima
		peak_means = argrelextrema(smooth_histo, np.greater)

		if not peak_means or len(peak_means) == 0:

			# Not enough data to create peaks in this case.
			self.peaks = None

		else:

			peaks = []

			# Unless otherwise specified we assume fmax represents the homozygous peak (primary content)
			# Explore expected peak sites, also look for heterozygous content.
			for mu in peak_means[0]:

				# Correct for smoothing
				mean = mu - wlen + 2

				# Take a first guess at the std dev, just use something relatively narrow for now
				# the local optimisation should pad this out later
				sigma = 2.0

				# We are (at present, only) interested in a region up to 2 stddevs from the mean (95% coverage)
				radius = int(sigma * 2.0)

				# Ignore anything too close to the edge
				if mean - radius > 0 and mean + radius < self.k:
					# This code assumes a maxima exists here
					peaks.append(Peak(
						mean,		# Mean
						sigma,  	# Std dev
						self.histogram[mean], # Use the histogram value at this position as an initial value for this peak
						mean == np.argmax(self.histogram)  # Whether or not this is the primary peak
					))

			self.peaks = peaks

		return