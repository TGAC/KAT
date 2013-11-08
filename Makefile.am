
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4

# DOXYGEN SUPPORT
include $(top_srcdir)/aminclude.am


# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = README.md COPYING AUTHORS NEWS doxygen.cfg

# FILES NOT TO INSTALL
dist_noinst_SCRIPTS = autogen.sh clean.sh

# SCRIPTS TO ADD TO BIN DIR
dist_bin_SCRIPTS = \
        scripts/kat_comp_reads.sh \
        scripts/dist_analysis.py

# SRC DIRS
SUBDIRS = src tests


