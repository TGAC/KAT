
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

SEQAN_ARCHIVE = deps/seqan-library-2.0.0.tar.gz

dist_pkgdata_DATA = $(SEQAN_ARCHIVE)

install-data-hook:
    cd $(DESTDIR)/deps; \
    cat $(SEQAN_ARCHIVE) | gunzip | $(TAR) -H gnu --list > uninstall_manifest.txt; \
    cat $(SEQAN_ARCHIVE) | gunzip | $(TAR) -H gnu --no-same-owner --extract; \
    rm --force $(SEQAN_ARCHIVE); \
    cat uninstall_manifest.txt | sed --expression='s/^\|$$/"/g' | xargs chmod a=rX,u+w

uninstall-local:
    cd $(DESTDIR)/deps; \
    cat uninstall_manifest.txt | sed --expression='s/ /\\ /g' | xargs rm --force; \
    rm --force uninstall_manifest.txt

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST =	README.md \
		COPYING \
		AUTHORS \
		INSTALL \
		ChangeLog \
		deps/seqan-library-2.0.0


# Scripts NOT to install
dist_noinst_SCRIPTS = autogen.sh clean.sh

# Scripts to install
dist_bin_SCRIPTS = scripts/dist_analysis.py

# SRC DIRS
SUBDIRS = deps/jellyfish-2.2.0 tests

AUTOMAKE_OPTIONS = subdir-objects

# Apply for all
CXXFLAGS = -g -O0 -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor -std=c++11 @AM_CXXFLAGS@
CPPFLAGS = -DCPLUSPLUS -I$(top_srcdir)/deps/seqan-library-2.0.0/include -I$(top_srcdir)/deps/jellyfish-2.2.0/include -I$(top_srcdir)/src/inc -I$(top_srcdir)/src @AM_CPPFLAGS@
LDFLAGS = @AM_LDFLAGS@
LIBS = $(top_srcdir)/deps/jellyfish-2.2.0/libjellyfish-2.0.la @AM_LIBS@

bin_PROGRAMS = kat

katdir = $(top_srcdir)
kat_HEADERS = \
		src/inc/matrix/matrix_metadata_extractor.hpp \
		src/jellyfish_helper.hpp \
		src/input_handler.hpp \
		src/plot_density.hpp \
		src/plot_spectra_cn.hpp \
		src/plot_spectra_mx.hpp \
		src/plot_spectra_hist.hpp \
		src/plot.hpp \
		src/comp.hpp \
		src/gcp.hpp \
		src/histogram.hpp \
		src/sect.hpp
		

kat_SOURCES = \
                src/inc/gnuplot/gnuplot_i.cc \
		src/inc/gnuplot/gnuplot_i.hpp \
		src/inc/matrix/sparse_matrix.hpp \
		src/inc/matrix/threaded_sparse_matrix.hpp \
                src/inc/matrix/matrix_metadata_extractor.cc \
                src/inc/str_utils.hpp \
		src/inc/spectra_helper.hpp \
                src/inc/kat_fs.hpp \
		src/jellyfish_helper.cc \
		src/input_handler.cc \
		src/plot_density.cc \
		src/plot_profile.hpp \
		src/plot_spectra_cn.cc \
		src/plot_spectra_mx.cc \
		src/plot_spectra_hist.cc \
		src/plot.cc \
		src/comp.cc \
		src/gcp.cc \
		src/histogram.cc \
		src/sect.cc \
		src/kat.cc


