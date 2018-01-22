AUTOMAKE_OPTIONS = subdir-objects


# What to build
bin_PROGRAMS = kat

# Executable inputs
kat_CXXFLAGS = \
	-g -O3 -fwrapv -Wall -Wextra \
	-Wno-deprecated-declarations -Wno-unused-function -Wno-unused-parameter \
	-Wno-unused-variable -Wno-unused-command-line-argument \
	@AM_CXXFLAGS@

kat_CPPFLAGS = \
	-isystem $(top_srcdir)/deps/seqan-library-2.0.0/include \
	-isystem $(top_srcdir)/deps/jellyfish-2.2.0/include \
	-isystem $(top_srcdir)/deps/boost_build/include \
	-isystem $(top_srcdir)/lib/include \
	-DKAT_SITE_PKGS='"$(libdir)/python$(PYTHON_VERSION)/site-packages"' \
	@AM_CPPFLAGS@

kat_LDFLAGS = \
	-static \
	-L../lib \
	-L../deps/jellyfish-2.2.0 \
	-L../deps/boost_build/lib \
	@AM_LDFLAGS@

kat_LDADD = \
	../deps/boost_build/lib/libboost_timer.a \
	../deps/boost_build/lib/libboost_chrono.a \
	../deps/boost_build/lib/libboost_filesystem.a \
	../deps/boost_build/lib/libboost_program_options.a \
	../deps/boost_build/lib/libboost_system.a \
	-lkat \
	-lkat_jellyfish \
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
