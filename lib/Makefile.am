AUTOMAKE_OPTIONS = subdir-objects

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = kat-2.1.pc

lib_LTLIBRARIES = libkat-2.1.la

libkat_2_1_la_LDFLAGS = -version-info 2:1:0
libkat_2_1_la_SOURCES = \
	src/gnuplot_i.cc \
	src/matrix_metadata_extractor.cc \
	src/input_handler.cc \
	src/jellyfish_helper.cc \
	src/comp_counters.cc

library_includedir=$(includedir)/kat-@PACKAGE_VERSION@/kat
KI = $(top_srcdir)/lib/include/kat
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

libkat_2_1_la_CPPFLAGS = \
	-isystem $(top_srcdir)/deps/jellyfish-2.2.0/include \
	-isystem $(top_srcdir)/lib/include \
	-DBINDIR=\"$(bindir)\" \
	@AM_CPPFLAGS@
	
libkat_2_1_la_CXXFLAGS = -g -O3 -fwrapv -Wall -Wextra -Wno-deprecated-declarations -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable -ansi -pedantic -Werror -std=c++11  @AM_CXXFLAGS@