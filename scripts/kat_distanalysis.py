#!/usr/bin/env python3

import argparse
from math import *
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
	def __init__(self, histo_file=None, points=10000, column=1, cumulative=False):
		"""Init with the objective histogram"""
		self.histogram = []
		self.peaks = []
		self.k = 27
		if histo_file:
			self.read_hist(histo_file, points, column, cumulative=cumulative)

	def read_hist(self, name, points=10000, column=1, cumulative=False):
		f = open(name)
		for l in f.readline():
			if l.startswith("# Kmer value:"):
				self.k = int(l.strip().split(":")[1])
				break

		if cumulative:
			self.histogram = [sum([int(y) for y in x.split()[column:]]) for x in f.readlines() if x and x[0] != '#'][
							 :points][1:]
		else:
			self.histogram = [int(x.split()[column]) for x in f.readlines() if x and x[0] != '#'][:points][1:]
		f.close()

	def total_values(self, start=1, end=10000):
		return list(map(sum, list(zip(*[x.points(start, end) for x in self.peaks]))))

	###----------------- FIND AND CREATE PEAKS --------------------

	def deriv(self, histogram):
		return [histogram[i + 2] - histogram[i] for i in range(len(histogram) - 2)]

	def progsmoothderiv(self, histogram):
		# always return mean
		return [mean(histogram[i + 1:int(i + 1 + i / 10)]) - mean(histogram[int(i - i / 10):i + 1]) for i in
				range(len(histogram) - 2)]

	def find_maxima(self, center, radius, min_perc, min_elem, histogram=None):
		if histogram == None:
			hs = self.histogram[center - radius:center + radius]
			# Smoothed df/dx (dx=2 units)
			deriv = self.progsmoothderiv(self.histogram)[center - radius:center + radius]
		else:
			hs = histogram[center - radius:center + radius]
			# Smoothed df/dx (dx=2 units)
			deriv = self.progsmoothderiv(histogram)[center - radius:center + radius]

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

	def create_peaks(self, min_perc, min_elem):
		fmin = 0
		# walk till first local minimum (d(f)>0)
		for i in range(len(self.histogram)):
			if self.histogram[i] < self.histogram[i + 1]:
				fmin = i
				break
		if not fmin:
			print(self.histogram[:25])
			print("Could not find local minima")
			self.fmax = 0
			self.maxval = 0

		# Find the global fm|p(fm)=max(p)
		fmax = self.histogram.index(max(self.histogram[fmin:]))
		# print "fmin=%d (%d) |fmax=%d (%d)" % (fmin,self.histogram[fmin],fmax,self.histogram[fmax])
		if fmax < 10:
			print("fmax<10, no analysis can be performed")
			return
		self.fmax = fmax
		self.fmin = fmin
		self.maxval = self.histogram[fmax]
		# Explore fm/x and fm*x for x in [1,4]
		for f in [fmax / 4, fmax / 3, fmax / 2, fmax, fmax * 2, fmax * 3, fmax * 4]:
			if f + f / 5 < len(self.histogram):
				lm = self.find_maxima(int(f), int(f / 5), min_perc, min_elem)
				if lm:
					self.add_peak_and_update_cuts(lm, reset_opt=True)
					# if not fdists:
					#    print "Local maxima not found, skipping further analysis"
					#    return
					# Guess counts for all relevant distributions
					# print "Local maxima found on %s" %(str(fdists))

	###----------------- OPTIMIZATION --------------------

	def optimize_peaks(self, min_perc, min_elem):
		print("------ New independent optimization started")
		sortedpeaks = [x for x in self.peaks]
		sortedpeaks.sort(key=lambda x: -x.elements)

		# create a base as the histogram and start from there
		base = [x for x in self.histogram]
		for p in sortedpeaks:
			i = self.peaks.index(p)
			# locally optimize the preak
			print("Optimizing peak: %s on [%d,%d]" % (p, self.cuts[i], self.cuts[i + 1]))
			p.constrained_opt(self.cuts[i], base[self.cuts[i]:self.cuts[i + 1]])
			print("   Result: %s" % p)
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
			print("Processing")
			min_peak_index = peak_index
			delta_peak_freq = 1000000
			for p_i, p in enumerate(self.peaks, start=1):
				delta = abs(p.mean - approx_freq)
				print("Delta", delta)
				if delta_peak_freq > delta:
					print("Mod")
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


class KmerSpectraAnalysis(object):
	def __init__(self, filename, points=10000):
		self.spectra = KmerSpectra(filename, points)

	def plot_hist(self, h, points, cap):
		plt.plot([min(cap, x) for x in h[:points]])

	def plot_hist_df(self, h, points, cap):
		plt.plot([max(-cap, min(cap, x)) for x in [h[i + 1] - h[i] for i in range(points)]])

	def plot_all(self, points, cap):
		self.plot_hist(self.spectra.histogram, points, cap)

		# self.plot_hist_df(self.histogram,points,cap)
		self.plot_hist(self.spectra.total_values(1, points + 1), points, cap)

	def plot_peaks(self, points, cap):
		for p in self.spectra.peaks:
			self.plot_hist(p.points(1, points + 1), points, cap)

	def analyse(self):
		self.spectra.create_peaks()
		limy = int(self.spectra.maxval * 1.1 / 1000) * 1000
		limx = self.spectra.peaks[-1].mean * 2
		print("Plot limits: y->%d, x->%d" % (limy, limx))
		self.plot_all(limx, limy)
		plt.show()
		self.spectra.optimize_peaks()
		plt.figure()
		self.plot_all(limx, limy)
		plt.show()
		self.spectra.optimize_overall()
		plt.figure()
		self.plot_all(limx, limy)
		plt.figure()
		self.plot_peaks(limx, limy)
		plt.show()
		for p in self.spectra.peaks: print(p)


class MXKmerSpectraAnalysis(object):
	def __init__(self, filename, columns=3, points=10000, hom_peak_freq=0):
		self.spectras = [KmerSpectra(filename, points, column=0, cumulative=True)]
		self.hom_peak = hom_peak_freq

		if self.hom_peak > 0:
			print("User specified homozygous peak frequency at:", self.hom_peak)
		for i in range(columns):
			self.spectras.append(KmerSpectra(filename, points, column=i, cumulative=(i == columns - 1)))

	def plot_hist(self, h, points, cap, label=""):
		plt.plot([min(cap, x) for x in h[:points]], label=label)
		plt.xlabel('Kmer Frequency')
		plt.ylabel('# Distinct Kmers')
		plt.legend()

	def plot_hist_df(self, h, points, cap):
		plt.plot([max(-cap, min(cap, x)) for x in [h[i + 1] - h[i] for i in range(points)]])

	def plot_all(self, points=0, cap=0, spectra=True, fit=True, dists=True):
		if 0 == points: points = self.limx
		if 0 == cap: cap = self.limy
		for s in self.spectras:
			if self.spectras[0] == s:
				slabel = "General Spectra"
			else:
				slabel = "%d x present" % (self.spectras.index(s) - 1)
			plt.figure()
			if spectra: self.plot_hist(s.histogram, points, cap, label=slabel)
			if fit: self.plot_hist(s.total_values(1, points + 1), points, cap, label=slabel + " fit")
			if dists:
				for p in s.peaks:
					self.plot_hist(p.points(1, points + 1), points, cap, label="fit dist %d" % s.peaks.index(p))
			plt.show()
			for p in s.peaks:
				print(p)

	def analyse(self, min_perc=1, min_elem=100000, verbose=False):
		self.limx = 0
		self.limy = 0
		for s in self.spectras:
			print("analysing spectra... ", end=' ')
			sys.stdout.flush()
			s.create_peaks(min_perc=min_perc, min_elem=min_elem)

			sys.stdout.flush()
			if s.peaks:
				self.limy = max(int(s.maxval * 1.1 / 1000) * 1000, self.limy)
				self.limx = max(min(s.peaks[-1].mean * 2, len(s.histogram)), self.limx)
				print("peaks created ... ", end=' ')
				sys.stdout.flush()
				s.optimize_peaks(min_perc=min_perc, min_elem=min_elem)
				print("locally optimised ... ", end=' ')
				for p in s.peaks: print(p)
				sys.stdout.flush()
				s.optimize_overall()
				print("overall optimised ... DONE")
				sys.stdout.flush()
		print("Plot limits: y->%d, x->%d" % (self.limy, self.limx))

	def peak_stats(self):
		"""TODO: Runs analyse (TODO:include general spectra)
				 Takes enough peaks as to cover a given % of the elements:
					 - Find the peak across all distributions
					 - Reports peak stats
				 If multiple peaks have been analyzed, tries to find the "main unique" and explains the results based on that freq.
		"""
		# step 1, try to find a reasonable mean for kmer frequency.
		# weighted means by number of elements?
		general_dists = [(x.mean, x.elements) for x in self.spectras[0].peaks]
		general_dists.sort(key=lambda x: x[0])
		pkc = general_dists[0][0]
		kcov = sum([(x[0] / round(x[0] / pkc)) * x[1] for x in general_dists]) / sum([x[1] for x in general_dists])
		print(pkc, kcov)
		# step 2, selects frequencies for peaks from bigger to smaller till X% of the elements are covered or no more peaks
		goal = 0.99 * sum([x[1] for x in general_dists])
		maxpeaks = 10
		general_dists.sort(key=lambda x: -x[1])
		af = []
		peaks = 0
		covered = 0
		for x in general_dists:
			af.append(kcov * round(x[0] / kcov))
			peaks += 1
			covered += x[1]
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
						" %dx: %.2f%% (%d elements on f=%.2f)" % (i, float(pd[i][1]) * 100 / total, pd[i][1], pd[i][0]))
				else:
					print(" %dx: No significant content" % i)

		# Step 4, genome stats
		hp = self.spectras[0].getHomozygousPeakIndex(self.hom_peak)
		gs = self.spectras[0].calcGenomeSize(hom_peak=hp)
		hr = self.spectras[0].calcHetRate(gs)
		print()
		print("K-value used:", self.spectras[0].k)
		print("Peaks in analysis:", len(self.spectras[0].peaks))
		print("Homozygous peak index:", hp)
		print("Estimated genome size:", '{0:.2f}'.format(float(gs) / 1000000.0), "Mbp")
		if (hp > 1):
			print("Heterozygous rate:", '{0:.2f}'.format(hr), "%")
		print()

		return




if __name__ == '__main__':
	# ----- command line parsing -----
	parser = argparse.ArgumentParser(
		description="""Analyse a comp matrix file with respect to the distributions and copy numbers seen within.""")

	parser.add_argument("matrix_file", type=str,
						help="The input matrix file from KAT comp")

	parser.add_argument("-c", "--cns", type=int, default=4,
						help="The number of copy numbers to consider in the analysis")
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


	args = parser.parse_args()

	a = MXKmerSpectraAnalysis(args.matrix_file, args.cns, args.freq_cutoff, args.homozygous_peak)
	a.analyse(min_perc=args.min_perc, min_elem=args.min_elem)
	a.peak_stats()
	if args.plot:
		a.plot_all()