
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# Blank man target (used in doc)
man: 
	cd doc && $(MAKE) $(AM_MAKEFLAGS) man
.PHONY: man
	
# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = \
	README.md \
	COPYING \
	AUTHORS \
	NEWS
	

# Scripts NOT to install
dist_noinst_SCRIPTS = autogen.sh antigen.sh

# Scripts to install
dist_bin_SCRIPTS = \
	scripts/kat_distanalysis.py \
	scripts/kat_plot_misc.py \
	scripts/kat_plot_colormaps.py \
	scripts/kat_plot_density.py \
	scripts/kat_plot_profile.py \
	scripts/kat_plot_spectra-cn.py \
	scripts/kat_plot_spectra-hist.py \
	scripts/kat_plot_spectra-mx.py

# SRC DIRS
make_dirs=deps/jellyfish-2.2.0 deps/seqan-library-2.0.0 lib src tests

if test @MAKE_DOCS@ = 1; then
    make_dirs +=  doc

SUBDIRS = $(make_dirs)

AUTOMAKE_OPTIONS = subdir-objects