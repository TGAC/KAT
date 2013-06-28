AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
EXTRA_DIST = README COPYING AUTHORS NEWS

AM_LDFLAGS = -pthread
AM_CPPFLAGS = -Wall -Wno-sign-compare -Wnon-virtual-dtor -I$(top_srcdir)/src/inc
AM_CXXFLAGS = -g -O3


bin_PROGRAMS = $(top_builddir)/bin/kat

__top_builddir__bin_kat_SOURCES = \
                            src/inc/gnuplot/gnuplot_i.cc \
                            src/inc/jellyfish/jellyfish_helper.hpp \
                            src/sect/sect_main.cc \
                            src/comp/comp_main.cc \
			    src/plot/plot_main.cc \
                            src/plot/flame/flame_plot_main.cc \
                            src/plot/asm/asm_plot_main.cc \
                            src/plot/sect/sect_plot_main.cc \
                            src/kat.cc

__top_builddir__bin_kat_LDADD = lib/libjellyfish-1.1.a lib/libz.a


dist_noinst_SCRIPTS = autogen.sh
dist_bin_SCRIPTS = scripts/kat_comp_reads.sh
