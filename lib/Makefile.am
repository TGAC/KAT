AUTOMAKE_OPTIONS = subdir-objects

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = kat.pc

lib_LTLIBRARIES = libkat.la

libkat_la_LDFLAGS = \
	-L$(top_builddir)/deps/boost/build/lib/ \
	-version-info 2:4:2

LIBS = \
	$(top_builddir)/deps/jellyfish-2.2.0/libkat_jellyfish.la \
	-lboost_timer \
	-lboost_chrono \
	-lboost_filesystem \
	-lboost_program_options \
	-lboost_system \
	-lz \
	@AM_LIBS@


libkat_la_SOURCES = \
	src/matrix_metadata_extractor.cc \
	src/input_handler.cc \
	src/jellyfish_helper.cc \
	src/comp_counters.cc

library_includedir=$(includedir)/kat-@PACKAGE_VERSION@/kat

KI = $(top_srcdir)/lib/include/kat
library_include_HEADERS =   $(KI)/distance_metrics.hpp \
			    $(KI)/input_handler.hpp \
			    $(KI)/jellyfish_helper.hpp \
			    $(KI)/kat_fs.hpp \
			    $(KI)/matrix_metadata_extractor.hpp \
			    $(KI)/sparse_matrix.hpp \
			    $(KI)/spectra_helper.hpp \
			    $(KI)/str_utils.hpp \
			    $(KI)/comp_counters.hpp

libkat_la_CPPFLAGS = \
	-I$(top_srcdir)/deps/jellyfish-2.2.0/include \
	-I$(top_srcdir)/deps/boost/build/include \
	-I$(top_srcdir)/lib/include \
	-DKAT_SITE_PKGS='"$(libdir)/python$(PYTHON_VERSION)/site-packages"' \
	@AM_CPPFLAGS@

libkat_la_CXXFLAGS = \
	-fwrapv -m64 \
	-Wall \
	-Wno-deprecated-declarations -Wno-unused-function -Wno-unused-parameter -Wno-unused-variable \
	@AM_CXXFLAGS@
