AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = kat

kat_noinstdir = $(top_srcdir)/src

kat_CXXFLAGS = -g -O3 -Wno-conversion-null -Wno-unused-result -Wno-unused-function -Wno-sign-compare -Wnon-virtual-dtor -std=c++11 @AM_CXXFLAGS@

kat_CPPFLAGS = -DCPLUSPLUS -I$(top_srcdir)/deps/seqan-library-2.0.0/include -I$(top_srcdir)/deps/jellyfish-2.2.0/include -I$(top_srcdir)/deps/jellyfish-2.2.0/include @AM_CPPFLAGS@

kat_LDADD = $(top_srcdir)/deps/jellyfish-2.2.0/.libs/libjellyfish-k2.0.a \
	    @AM_LIBS@

noinst_HEADERS = \
		inc/matrix/matrix_metadata_extractor.hpp \
		jellyfish_helper.hpp \
		input_handler.hpp \
		plot_density.hpp \
		plot_spectra_cn.hpp \
		plot_spectra_mx.hpp \
		plot_spectra_hist.hpp \
		plot.hpp \
		comp.hpp \
		gcp.hpp \
		histogram.hpp \
		sect.hpp
		

kat_SOURCES = \
		inc/gnuplot/gnuplot_i.cc \
		inc/gnuplot/gnuplot_i.hpp \
		inc/matrix/sparse_matrix.hpp \
		inc/matrix/threaded_sparse_matrix.hpp \
                inc/matrix/matrix_metadata_extractor.cc \
                inc/str_utils.hpp \
		inc/spectra_helper.hpp \
                inc/kat_fs.hpp \
		jellyfish_helper.cc \
		input_handler.cc \
		plot_density.cc \
		plot_profile.hpp \
		plot_spectra_cn.cc \
		plot_spectra_mx.cc \
		plot_spectra_hist.cc \
		plot.cc \
		comp.cc \
		gcp.cc \
		histogram.cc \
		sect.cc \
		kat.cc


