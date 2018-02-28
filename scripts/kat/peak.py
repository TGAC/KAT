import numpy as np
from scipy import optimize


def gaussian(x, mu, sig):
	return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def createModel(X, mu, sig, scale):
	model = np.zeros_like(X)
	for i, x in enumerate(X):
		model[i] = gaussian(x, mu, sig) * scale
	return model


class Peak(object):
	"""
	A distribution representing kmers covered a certain number of times.
	Contains methods for fitting to an interval
	"""

	def __init__(self, mean, stddev, peak, primary, description=""):
		self._mean = mean
		self._stddev = stddev
		self._peak = peak
		self.primary = primary
		self.Tx = None
		self.Ty = None
		self.description = description

	def left(self):
		return self._mean - self.radius()

	def right(self):
		return self._mean + self.radius()

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
		return self._stddev

	def mean(self, mean=None):
		if mean is not None:
			self._mean = mean
		return self._mean

	def peak(self, peak=None):
		if peak is not None:
			self._peak = peak
		return self._peak

	def elements(self):
		return int(sum(self.Ty)) if self.Ty is not None else 0

	# def poisson(self, x):
	#    return np.exp(-self._mean) * (np.power(self._mean, x) / math.factorial(x))

	def gaussian(self, x):
		return np.exp(-np.power(x - self._mean, 2.) / (2 * np.power(self._stddev, 2.)))

	def __str__(self):
		return "Peak of " + str(int(self._peak)) + " at frequency " + "{:.2f}".format(self._mean) + "(stddev: " + "{:.2f}".format(self._stddev) + "), with volume of " + \
			   str(self.elements()) + " elements between frequencies of " + "{:.2f}".format(self.left()) + " and " + "{:.2f}".format(
			self.right()) + "; Primary: " + str(self.primary)

	def toRow(self):
		return ["{:.2f}".format(self.left()), "{:.2f}".format(self._mean), "{:.2f}".format(self.right()), "{:.2f}".format(self._stddev),
			 str(int(self._peak)), str(int(self.elements())), str(self.description)]

	@staticmethod
	def header():
		return ["Left","Mean","Right","StdDev","Max","Volume","Description"]

	def updateModel(self, new_mean, new_peak, new_stddev):
		"""
		Updates the histogram representing the gaussian modelled by this peak
		"""
		self._mean = new_mean
		self._peak = new_peak
		self._stddev = new_stddev

		for i, x in enumerate(self.Tx):
			self.Ty[i] = gaussian(x, self._mean, self._stddev) * self._peak

		return self.Ty

	def residuals(self, p, fmin=0):
		"""
		Fit this gaussian distribution as closely as possible to the histogram using the given parameters
		:param p: The parameters to use
		:return: Value representing the delta between the fitted gaussian and the histogram
		"""

		# This set the peak and adjusts the scaling factor accordingly
		# and then updates the histogram represented by this specific peak
		model = np.zeros_like(self.Tx)
		for i, x in enumerate(self.Tx):
			model[i] = gaussian(x, p[0], p[2]) * p[1]

		# Return the distance between the fitted peak and the actual histogram at each site
		residuals = self.histogram - model

		# We want to more heavily penalise all points which exceed the histogram
		for i in range(len(residuals)):
			d = residuals[i]
			#if d < 0:
			#	residuals[i] = d * 100

			# Suppress residuals that come before fmin (we aren't interested in fitting to the error K=mers)
			if i <= fmin:
				residuals[i] /= np.power(fmin - i + 1, 10)

		# The residual differences between the actual histogram and the model peak at each value of X
		return residuals

	def optimise(self, histogram, fmin=0):
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

		# Make sure the current settings are up to date and the fitted histogram is based on these
		self.updateModel(self._mean, self._peak, self._stddev)

		# Set up variables to optimise (just the peak)
		p = []
		lower_bounds = []
		upper_bounds = []
		p.append(self._mean)
		lower_bounds.append(self._mean - 1.0)
		upper_bounds.append(self._mean + 1.0)
		p.append(self._peak.astype(np.float64))
		lower_bounds.append(0.0)
		upper_bounds.append(self._peak.astype(np.float64))
		p.append(self._stddev)
		lower_bounds.append(1.0)
		upper_bounds.append(max((self._mean - 2.0) / 2.0, self._stddev))

		# Set the optimal peak value that maximises the space under the histogram, without going over the borders.
		res = optimize.least_squares(self.residuals, np.array(p).astype(np.float64), args=[fmin], bounds=(lower_bounds, upper_bounds), loss="soft_l1")

		# If all went well update the model with the optimised variables
		if res.success:
			self.updateModel(res.x[0], res.x[1], res.x[2])
		else:
			raise ValueError("Problem optimising peak.")

		return
