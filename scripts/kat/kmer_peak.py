import numpy as np
import math
from scipy import optimize


class KmerPeak(object):
	"""
	A distribution representing kmers covered a certain number of times.
	Contains methods for fitting to an interval
	"""

	def __init__(self, mean, stddev, peak, primary):
		self._mean = mean
		self._stddev = stddev
		self._peak = peak
		self.primary = primary
		self._scaling_factor = self.calcScalingFactor()
		self.Tx = None
		self.Ty = None

	def left(self):
		return int(self._mean - self.radius())

	def right(self):
		return int(self._mean + self.radius())

	def radius(self):
		"""
		Returns the radius we are interested in for this peak.
		Current set at 2 * the stddev, so 2 * radius gives 95% coverage
		:return: region of interest from mean
		"""
		return 2.0 * self._stddev

	def stddev(self, stddev=None):
		if stddev is not None:
			self._stddev = stddev
			self._scaling_factor = self.calcScalingFactor()
		return self._stddev

	def mean(self, mean=None):
		if mean is not None:
			self._mean = mean
			self._scaling_factor = self.calcScalingFactor()
		return self._mean

	def peak(self, peak=None):
		if peak is not None:
			self._peak = peak
			self._scaling_factor = self.calcScalingFactor()
		return self._peak

	def elements(self):
		return int(sum(self.Ty)) if self.Ty is not None else 0

	# def poisson(self, x):
	#    return np.exp(-self._mean) * (np.power(self._mean, x) / math.factorial(x))

	def gaussian(self, x):
		return np.exp(-np.power(x - self._mean, 2.) / (2 * np.power(self._stddev, 2.)))

	def calcScalingFactor(self):
		p = self.gaussian(self._mean)
		return (float(self._peak) / p) if p > 0 else 0

	def __str__(self):
		return "Peak of " + str(self._peak) + " at frequency " + str(self._mean) + ", with volume of " + \
			   str(self.elements()) + " elements between frequencies of " + str(self.left()) + " and " + str(
			self.right()) + "; Primary: " + str(self.primary)

	def toTabString(self):
		return "\t".join(
			[str(self.left()), str(int(self._mean)), str(self.right()), str(int(self._peak)), str(int(self.elements())),
			 str(self.primary)])

	@staticmethod
	def getTabHeader():
		return "Left\tMean\tRight\tMax\tVolume\tPrimary"

	def updateModel(self, new_peak, new_stddev):
		"""
		Updates the histogram representing the gaussian modelled by this peak
		"""
		self._peak = new_peak
		self._stddev = new_stddev
		self._scaling_factor = self.calcScalingFactor()

		# TODO there's probably a super fast numpy vectorised way of doing this.
		for i, x in enumerate(self.Tx):
			self.Ty[i] = int(self.gaussian(x) * self._scaling_factor)

		return self.Ty

	def _objectiveFunc(self, p):
		"""
		Fit this gaussian distribution as closely as possible to the histogram using the given parameters
		:param p: The parameters to use
		:return: Value representing the delta between the fitted gaussian and the histogram
		"""

		# This set the peak and adjusts the scaling factor accordingly
		# and then updates the histogram represented by this specific peak
		self.updateModel(p[0], p[1])

		# Return the distance between the fitted peak and the actual histogram at each site
		delta = np.abs(self.Ty - self.histogram)

		# We want to heavily penalise all points which exceed the histogram
		#for i in range(len(delta)):
		#	d = delta[i]
		#	delta[i] = np.power(d + 1000, 2) if d > 0 else -d

		# Sum the distances to provide overall level of difference
		return sum(delta)

	def optimise(self, histogram):
		"""
		Tries to fit this single guassian distribution to this point in the histogram as closely as possible
		:param histogram:
		:return:
		"""

		# Sanity check...
		if len(histogram) == 0:
			raise RuntimeError("Can't model")

		# Save histogram for easy access
		self.histogram = np.array(histogram)

		# Create Tx and Ty (represents fitted histograms)
		self.Tx = np.linspace(0, len(histogram) - 1, len(histogram))
		self.Ty = np.zeros_like(self.Tx)

		# Update the fitted histograms based on this peak's scaled gaussian
		self.updateModel(self._peak, self._stddev)

		# Set up variables to optimise (just the peak)
		p = [float(self._peak), self._stddev]
		bounds = [(0.1, self.histogram[self._mean]), (np.sqrt(np.sqrt(self._mean)), self._mean)]

		# Set the optimal peak value that maximises the space under the histogram, without going over the borders.
		res = optimize.minimize(self._objectiveFunc, p, bounds=bounds)
		if not res.success:
			print(
				"It is likely that the spectra is too complex to analyse properly.  Stopping analysis.\nOptimisation results:\n" + str(
					res[-2]))

		return
