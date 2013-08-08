AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}


# additional include paths necessary to compile the C++ library
AM_CXXFLAGS = -I$(top_srcdir)/src @AM_CXXFLAGS@

AM_LDFLAGS = -all-static -pthread
AM_CPPFLAGS = -Wall -Wno-sign-compare -Wnon-virtual-dtor -Iinc

bin_PROGRAMS = kat

kat_SOURCES = \
                inc/gnuplot/gnuplot_i.cc \
                inc/matrix/matrix_metadata_extractor.cc \
                sect/sect_main.cc \
                comp/comp_main.cc \
                gcp/gcp_main.cc \
                plot/plot_main.cc \
                plot/flame/flame_plot_main.cc \
                plot/asm/asm_plot_main.cc \
                plot/sect/sect_plot_main.cc \
                kat.cc

