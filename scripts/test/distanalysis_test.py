import unittest
import os
import tempfile
import shutil

from kat.distanalysis import *

class DistAnalysisTest(unittest.TestCase):

	@classmethod
	def setUpClass(self):
		self.res_dir = os.path.join(os.path.dirname(__file__), "resources")
		self.temp_dir = tempfile.mkdtemp(dir=self.res_dir)

	@classmethod
	def tearDownClass(self):
		shutil.rmtree(self.temp_dir, ignore_errors=True)

	def test_system_hist1(self):
		in_file = os.path.join(os.path.dirname(__file__), "resources", "hist1.hist")
		a = HistKmerSpectraAnalysis(in_file, haploid=False, freq_cutoff=500, k=27)
		a.analyse()
		a.peak_stats(os.path.join(self.temp_dir, "system_hist1"))
		assert(os.path.exists(os.path.join(self.temp_dir, "system_hist1.dist_analysis.json")))


	def test_system_gcp1(self):
		in_file = os.path.join(os.path.dirname(__file__), "resources", "gcp1.mx")
		a = GCKmerSpectraAnalysis(in_file, haploid=False, freq_cutoff=500, k=27)
		a.analyse()
		a.peak_stats(os.path.join(self.temp_dir, "system_gcp1"))
		assert(os.path.exists(os.path.join(self.temp_dir, "system_gcp1.dist_analysis.json")))


	def test_system_spectracn1(self):
		in_file = os.path.join(os.path.dirname(__file__), "resources", "spectracn1.mx")
		a = MXKmerSpectraAnalysis(in_file, haploid=False, freq_cutoff=500, k=27)
		a.analyse()
		a.peak_stats(os.path.join(self.temp_dir, "system_spectracn1"))
		assert (os.path.exists(os.path.join(self.temp_dir, "system_spectracn1.dist_analysis.json")))


	def test_system_spectracn2(self):
		in_file = os.path.join(os.path.dirname(__file__), "resources", "spectracn2.mx")
		a = MXKmerSpectraAnalysis(in_file, haploid=False, freq_cutoff=500, k=27)
		a.analyse()
		a.peak_stats(os.path.join(self.temp_dir, "system_spectracn2"))
		assert (os.path.exists(os.path.join(self.temp_dir, "system_spectracn2.dist_analysis.json")))
