
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = kat-2.1.pc

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST =	README.md \
		COPYING \
		AUTHORS \
		INSTALL \
		ChangeLog \
		tests/data/EcoliK12.fasta \
		tests/data/ecoli.header.jf27 \
		tests/data/kat.hist \
		tests/data/sect_length_test.fa \
		tests/data/sect_test.fa \
		tests/data/ecoli_r1.1K.fastq \
		tests/data/ecoli_r2.1K.fastq \
		tests/test_comp.sh \
		tests/test_gcp.sh \
		tests/test_hist.sh \
		tests/test_sect.sh
	

# Scripts NOT to install
dist_noinst_SCRIPTS = autogen.sh antigen.sh

# Scripts to install
dist_bin_SCRIPTS =	scripts/kat_distanalysis.py \
			scripts/kat_plot_misc.py \
			scripts/kat_plot_colormaps.py \
			scripts/kat_plot_density.py \
			scripts/kat_plot_profile.py \
			scripts/kat_plot_spectra-cn.py \
			scripts/kat_plot_spectra-hist.py \
			scripts/kat_plot_spectra-mx.py

# What to build
bin_PROGRAMS = src/kat
lib_LTLIBRARIES = lib/libkat-2.1.la


# Executable inputs
src_kat_CXXFLAGS = -g -O3 -fwrapv -Wall -Wextra -Wno-deprecated-declarations -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -ansi -pedantic -Werror -std=c++11 @AM_CXXFLAGS@
src_kat_CPPFLAGS = -isystem $(srcdir)/deps/seqan-library-2.0.0/include -isystem $(srcdir)/deps/jellyfish-2.2.0/include -isystem $(srcdir)/include @AM_CPPFLAGS@

src_kat_LDADD = deps/jellyfish-2.2.0/.libs/libkat_jellyfish-2.0.la \
		lib/libkat-2.1.la \
		@AM_LIBS@

noinst_HEADERS = \
		src/plot_density.hpp \
		src/plot_spectra_cn.hpp \
		src/plot_spectra_mx.hpp \
		src/plot_spectra_hist.hpp \
		src/plot.hpp \
		src/filter_kmer.hpp \
		src/filter_sequence.hpp \
		src/filter.hpp \
		src/comp.hpp \
		src/gcp.hpp \
		src/histogram.hpp \
		src/sect.hpp
		

src_kat_SOURCES = \
		src/plot_density.cc \
		src/plot_profile.hpp \
		src/plot_spectra_cn.cc \
		src/plot_spectra_mx.cc \
		src/plot_spectra_hist.cc \
		src/plot.cc \
		src/filter_kmer.cc \
		src/filter_sequence.cc \
		src/filter.cc \
		src/comp.cc \
		src/gcp.cc \
		src/histogram.cc \
		src/sect.cc \
		src/kat.cc


# Library inputs
lib_libkat_2_1_la_LDFLAGS = -version-info 2:1:0
lib_libkat_2_1_la_SOURCES = lib/gnuplot_i.cc \
			lib/matrix_metadata_extractor.cc \
			lib/input_handler.cc \
			lib/jellyfish_helper.cc \
			lib/comp_counters.cc

library_includedir=$(includedir)/kat-@PACKAGE_VERSION@/kat
KI = include/kat
library_include_HEADERS =   $(KI)/distance_metrics.hpp \
			    $(KI)/gnuplot_i.hpp \
			    $(KI)/input_handler.hpp \
			    $(KI)/jellyfish_helper.hpp \
			    $(KI)/kat_fs.hpp \
			    $(KI)/matrix_metadata_extractor.hpp \
			    $(KI)/sparse_matrix.hpp \
			    $(KI)/spectra_helper.hpp \
			    $(KI)/str_utils.hpp \
			    $(KI)/comp_counters.hpp

lib_libkat_2_1_la_CPPFLAGS = -isystem $(srcdir)/deps/jellyfish-2.2.0/include -isystem $(srcdir)/include @AM_CPPFLAGS@


# SRC DIRS
SUBDIRS = deps/jellyfish-2.2.0 deps/seqan-library-2.0.0 tests

AUTOMAKE_OPTIONS = subdir-objects