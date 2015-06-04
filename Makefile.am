
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = include \
README.md \
COPYING \
AUTHORS \
deps/jellyfish-2.2.0/**.cc \
deps/jellyfish-2.2.0/**.cpp \
deps/jellyfish-2.2.0/**.hpp \
deps/jellyfish-2.2.0/**.h \
deps/jellyfish-2.2.0/**.i \
deps/jellyfish-2.2.0/**.py \
deps/jellyfish-2.2.0/**.pm \
deps/jellyfish-2.2.0/**.t \
deps/jellyfish-2.2.0/**.am \
deps/jellyfish-2.2.0/**.sh \
deps/jellyfish-2.2.0/**.in \
deps/jellyfish-2.2.0/**.m4 \
deps/jellyfish-2.2.0/doc/* \
deps/jellyfish-2.2.0/README \
deps/jellyfish-2.2.0/LICENSE \
deps/jellyfish-2.2.0/configure \
deps/jellyfish-2.2.0/config.guess \
deps/jellyfish-2.2.0/configure.ac \
deps/jellyfish-2.2.0/depcomp \
deps/jellyfish-2.2.0/gtest.mk \
deps/jellyfish-2.2.0/install-sh \
deps/jellyfish-2.2.0/missing \
deps/jellyfish-2.2.0/test-driver

# FILES NOT TO INSTALL
dist_noinst_SCRIPTS = autogen.sh clean.sh

# SRC DIRS
SUBDIRS = deps/jellyfish-2.2.0

AUTOMAKE_OPTIONS = subdir-objects

# Apply for all
CXXFLAGS = -g -O0 -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor -std=c++11 @AM_CXXFLAGS@
CPPFLAGS = -DCPLUSPLUS -I$(top_srcdir)/deps/seqan-library-2.0.0/include -I$(top_srcdir)/deps/jellyfish-2.2.0/include -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @AM_CPPFLAGS@
LDFLAGS = @AM_LDFLAGS@
LIBS = $(top_srcdir)/deps/jellyfish-2.2.0/libjellyfish-2.0.la @AM_LIBS@

bin_PROGRAMS = kat

kat_SOURCES = \
                src/inc/gnuplot/gnuplot_i.cc \
		src/inc/gnuplot/gnuplot_i.hpp \
		src/inc/fstream_default.hpp \
		src/inc/matrix/sparse_matrix.hpp \
		src/inc/matrix/threaded_sparse_matrix.hpp \
		src/inc/matrix/matrix_metadata_extractor.hpp \
                src/inc/matrix/matrix_metadata_extractor.cc \
                src/inc/str_utils.hpp \
		src/inc/spectra_helper.hpp \
                src/inc/kat_fs.hpp \
		src/inc/kat_hash_counter.hpp \
                src/jellyfish_helper.cc \
		src/input_handler.cc \
		src/comp.cc \
		src/gcp.cc \
		src/histogram.cc \
		src/plot.hpp \
		src/plot_density.hpp \
		src/plot_spectra_cn.hpp \
                src/plot_profile.hpp \
		src/plot_spectra_hist.hpp \
		src/plot_spectra_mx.hpp \
                src/sect.cc \
		src/kat.cc


