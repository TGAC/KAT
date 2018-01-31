AUTOMAKE_OPTIONS = subdir-objects


# What to build
bin_PROGRAMS = kat

# Executable inputs
kat_CXXFLAGS = \
	-g -O3 \
	-fwrapv -Wall -Wextra \
	-Wno-deprecated-declarations -Wno-unused-function -Wno-unused-parameter \
	-Wno-unused-variable -Wno-unused-command-line-argument \
	@AM_CXXFLAGS@

kat_CPPFLAGS = \
	-I$(top_srcdir)/deps/seqan-library-2.0.0/include \
	-I$(top_srcdir)/deps/jellyfish-2.2.0/include \
	-I$(top_srcdir)/deps/boost/build/include \
	-I$(top_srcdir)/lib/include \
	-DKAT_SITE_PKGS='"$(libdir)/python$(PYTHON_VERSION)/site-packages"' \
	@AM_CPPFLAGS@

kat_LDFLAGS = \
	-static \
	-L../lib \
	-L../deps/jellyfish-2.2.0 \
	-L../deps/boost/build/lib \
	@AM_LDFLAGS@

kat_LDADD = \
	-lkat \
	-lkat_jellyfish \
	../deps/boost/build/lib/libboost_timer.a \
	../deps/boost/build/lib/libboost_chrono.a \
	../deps/boost/build/lib/libboost_filesystem.a \
	../deps/boost/build/lib/libboost_program_options.a \
	../deps/boost/build/lib/libboost_system.a \
	@AM_LIBS@

noinst_HEADERS = \
	plot.hpp \
	filter_kmer.hpp \
	filter_sequence.hpp \
	filter.hpp \
	comp.hpp \
	gcp.hpp \
	histogram.hpp \
	sect.hpp

kat_SOURCES = \
	plot.cc \
	filter_kmer.cc \
	filter_sequence.cc \
	filter.cc \
	comp.cc \
	gcp.cc \
	histogram.cc \
	sect.cc \
	kat.cc
