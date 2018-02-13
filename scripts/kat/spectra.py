import abc
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.signal import argrelextrema
import tabulate

try:
	from peak import Peak, gaussian, createModel
except:
	from kat.peak import Peak, gaussian, createModel


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
	def _createInitialPeaks(self):
		pass

	def _updateModel(self, params):
		"""
		This function updates the fitted histogram based on the current parameters in each of the peaks in this spectra
		:return: The newly fitted histogram (self.fitted_histogram)
		"""

		if len(params) != len(self.peaks) * 3:
			raise ValueError("Unexpected number of parameters")

		for i in range(len(self.peaks)):
			new_mean = params[i * 3]
			new_peak = params[i * 3 + 1]
			new_stddev = params[i * 3 + 2]
			self.peaks[i].updateModel(new_mean, new_peak, new_stddev)

		self.Ty = np.zeros_like(self.Tx)
		for p in self.peaks:
			self.Ty += p.Ty

		return self.Ty

	def _createModel(self, x, *params):
		"""
		This creates a model based on the parameters provided.  We expect these to be a multiple of the number of peaks.
		The parameters control how the shape of the scaled guassians representing each peak.
		:params x: The x values, each of these will be applied to this function to create y which is of the same length
		:param params: New set of parameters adjusted by the optimiser
		:return: A numpy array of scalar values representing the difference between the model and reality at each X value
		"""

		if len(params) != len(self.peaks) * 3:
			raise ValueError("Unexpected number of parameters")

		# Recalculate the model based on information in the new parameters
		y = np.zeros_like(x)
		for i, peak in enumerate(self.peaks):
			new_mean = params[i * 3]
			new_peak = params[i * 3 + 1]
			new_stddev = params[i * 3 + 2]
			pdist = createModel(x, new_mean, new_stddev, new_peak)
			y += pdist

		return y

	def optimise(self, fmin=0):
		"""
		Given the full set of peaks, adjust all their heights in order to best fit the acutal histogram
		We also put some bounds around the limits these values can take to stop them going crazy
		"""
		if not self.peaks:
			print("Can't optimise peaks because none are defined.", end="", flush=True)
			return

		params = []
		lower_bounds = []
		upper_bounds = []
		for p in self.peaks:
			params.append(p.mean())
			lower_bounds.append(p.mean() - 2.0)
			upper_bounds.append(p.mean() + 2.0)
			params.append(p.peak())
			lower_bounds.append(0.0)
			upper_bounds.append(p.peak())
			params.append(p.stddev())
			stddev_lower = p.stddev() - np.sqrt(p.stddev())
			stddev_upper = max(min((p.mean() - 2.0) / 2.0, p.stddev() + np.sqrt(p.stddev())), p.stddev()+0.01)
			lower_bounds.append(stddev_lower)	# Make sure we can't get massively smaller
			upper_bounds.append(stddev_upper)	# Make sure we can't get too much bigger or make peak extend past 0 freq

		# Reset Tx and Ty in case the histogram has been modified
		self.Tx = np.linspace(0, len(self.histogram) - 1, len(self.histogram))

		# Suppress error k-mers
		fitcurve = np.array(self.histogram)
		for i in range(len(fitcurve)):
			if i <= fmin:
				fitcurve[i] /= np.power(fmin - i + 1, 6)

		# Fit model to real histogram
		res = optimize.curve_fit(self._createModel, self.Tx, fitcurve, p0=params, bounds=(np.array(lower_bounds), np.array(upper_bounds)))

		# Update the model with the optimised variables
		self._updateModel(res[0])
		return


	def analyse(self, min_elements=1, verbose=False):
		"""
		Analyse the histogram for peaks
		:param verbose: Prints additional information about progress to stdout
		"""

		if verbose:
			print()
			print("Creating initial peaks ... ", end="", flush=True)
		self._createInitialPeaks()

		if self.peaks:
			if verbose:
				print("done.", len(self.peaks), "peaks initially created")
				print()
				self.printPeaks()
				print()
				print("Locally optimising each peak ... ", end="")
			for p_i, p in enumerate(self.peaks):
				p.optimise(self.histogram)

			# For debugging
			if False:
				plt.plot(self.histogram, color='black')
				for p_i, p in enumerate(self.peaks):
					plt.plot(p.Ty)
				plt.xlim(0, 70)
				plt.ylim(0, 200000000)
				plt.show()

			# Remove any peaks that contain little to no content
			self.peaks = list(filter(lambda p: p.elements() >= min_elements, self.peaks))

			if verbose:
				print("done.")
				print()
				self.printPeaks()
				print()
				print("Fitting cumulative distribution to histogram by adjusting peaks ... ", end="", flush=True)
			try:
				self.optimise(fmin=self.fmin if type(self) == KmerSpectra else 0)
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
			header = ["Index"] + Peak.header()
			rows = [[str(p_i)] + p.toRow() for p_i, p in enumerate(self.peaks, start=1)]
			print(tabulate.tabulate(rows, header))
		else:
			print("No peaks detected")

	def plot(self, xmax, ymax, title=None, to_screen=True, output_file=None):

		fig = plt.figure()

		labels = []
		labels.append(plt.plot(self.histogram[:xmax], label="Actual", color="black"))
		for p in self.peaks:
			colour=None
			if p.description.startswith("1X"):
				colour="red"
			elif p.description.startswith("1/2X"):
				colour="blue"
			elif p.description.startswith("2X"):
				colour="green"
			elif p.description.startswith("3X"):
				colour = "orange"
			labels.append(plt.plot(p.Ty[:xmax], label=p.description, color=colour))
		labels.append(plt.plot(self.Ty[:xmax], label="Fitted model", color="gray"))

		plt.xlabel('Kmer Frequency' if type(self) == KmerSpectra else 'GC count')
		plt.ylabel('# Distinct Kmers')
		if title:
			plt.title(title)
		plt.xlim((0, xmax))
		plt.ylim((0, ymax))
		plt.legend()

		if to_screen:
			plt.show()

		if output_file:
			fig.savefig(output_file)






class KmerSpectra(Spectra):
	"""
	A kmer spectra, comprised of different peaks. Contains the general fitting method.
	"""

	def __init__(self, histogram, haploid=False, k=27):
		"""
		Inititalise the spectra with the actual histogram to model
		:param histogram: Histogram derived from one of the KAT tools
		:param k: K value used to construct the histogram
		"""

		# Initialise super
		Spectra.__init__(self, histogram, k)

		# Extra properties for K-mer spectra
		self.haploid = haploid	# If haploid then we don't look for the heterozygous (1/2) peak
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
			frequencies = []
			if not self.haploid:
				frequencies.append(fmax / 2.0)
			for i in range(1, 5):
				frequencies.append(fmax * i)

			for mu in frequencies:

				# In a poisson distribution the mean is the variance, so the stddev is simply the square root of the mean
				sigma = np.sqrt(mu)

				# We are (at present, only) interested in a region up to 2 stddevs from the mean (95% coverage)
				radius = int(sigma * 2.0)
				mean = int(mu)

				# Conditions:
				# - we need at least a radius of 2
				# - the peak frequency must be greater than fmin
				# - The extent of the distribution including the radius should extend over the histogram limits
				# - the peak size must be at least 1
				if radius >= 2 and mean > fmin and mu - radius > 0 and mu + radius < len(self.histogram) and self.histogram[mean] >= 1:
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
				if abs(p.mean() - self.fmax) < 1.0:
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
			if p_i >= hom_peak_index:
				delta = p_i - hom_peak_index + 1
				sum += delta * p.elements()
				p.description = str(delta) + "X"
			elif p_i < hom_peak_index:
				delta = hom_peak_index - p_i + 1
				sum += p.elements() / delta
				p.description = "1/" + str(delta) + "X"

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

		print()
		self.printPeaks()
		print()
		print("K-value used:", self.k)
		print("Peaks in analysis:", len(self.peaks))
		print("Global minima @ Frequency=" + str(self.fmin) + " (" + str(self.histogram[self.fmin]) + ")")
		print("Global maxima @ Frequency=" + str(self.fmax) + " (" + str(self.histogram[self.fmax]) + ")")
		print("Overall mean k-mer frequency:", str(self.calcKmerCoverage()) + "x")
		print()
		print("Calculating genome statistics")
		print("-----------------------------")
		if hom_peak_freq > 0:
			print("User-specified that homozygous peak should have a frequency of", hom_peak_freq)
		else:
			print("Assuming that homozygous peak is the largest in the spectra with frequency of:", int(self.peaks[hp].mean()))
		print("Homozygous peak index:", hp)
		print("CAUTION: the following estimates are based on having a clean spectra and having identified the correct homozygous peak!")
		print("Estimated genome size:", '{0:.2f}'.format(float(gs) / 1000000.0), "Mbp")
		if (hp > 1):
			hr = self.calcHetRate(gs)
			print("Estimated heterozygous rate:", "{0:.2f}".format(hr) + "%")


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