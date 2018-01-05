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
    version="2.4.0",
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
    entry_points={"console_scripts": ["kat_distanalysis = kat.kat_distanalysis:main",
                                      "kat_plot_density = kat.kat_plot_density:main",
                                      "kat_plot_profile = kat.kat_plot_profile:main",
                                      "kat_plot_spectra_cn = kat.kat_plot_spectra_cn:main",
                                      "kat_plot_spectra_hist = kat.kat_plot_spectra_hist:main",
                                      "kat_plot_spectra_mx = kat.kat_plot_spectra_mx:main"]},
    install_requires=['numpy','scipy','matplotlib']
)
