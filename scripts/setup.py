# coding: utf-8

"""Setup file for PyPI"""

from setuptools import setup, find_packages
from codecs import open
from os import path
import glob
import re
import sys

if sys.version_info.major != 3:
    raise EnvironmentError("""kat is a python module that requires python3,
    and is not compatible with python2.""")

setup(
    name="kat",
    version="2.4.1",
    description="A set of python scripts to support the Kmer Analysis Toolkit",
    long_description='''The Kmer Analysis Toolkit is a suite of tools that analyse
                        jellyfish hashes or sequence files (fasta or fastq) using kmer counts.
                        This python package provides additional support and utilities for
                        that program.  In particular, plotting functionality.''',
    url="https://github.com/TGAC/KAT",
    author="Daniel Mapleson",
    author_email="daniel.mapleson@earlham.ac.uk",
    license="GPLV3",
    zip_safe=False,
    keywords="kmer genomics",
    packages=find_packages(),
    entry_points={"console_scripts": ["kat_distanalysis = kat.distanalysis:main",
                                      "kat_plot_density = kat.plot.density:main",
                                      "kat_plot_profile = kat.plot.profile:main",
                                      "kat_plot_spectra_cn = kat.plot.spectra_cn:main",
                                      "kat_plot_spectra_hist = kat.plot.spectra_hist:main",
                                      "kat_plot_spectra_mx = kat.plot.spectra_mx:main",
                                      "kat_plot_cold = kat.plot.cold:main"]},
    test_suite="nose.collector",
	tests_require = [
		'nose',
    ],
	install_requires=['numpy','scipy','matplotlib','tabulate']
)
