AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

AM_LDFLAGS = -all-static -pthread
AM_CPPFLAGS = -Wall -Wno-sign-compare -Wnon-virtual-dtor -Iinc
AM_CXXFLAGS = -g -O3

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

