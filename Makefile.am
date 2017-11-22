
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

#BOOST
foo.c: foo.in
        FOOPLACE=`./thescript`; $(FOOPLACE) -c -o $@ $<

#git submodule update --init tools/build/ tools/auto_index/ tools/bcp tools/boostdep tools/inspect tools/litre
#git submodule update --init libs/program_options libs/chrono libs/timer libs/filesystem libs/system libs/stacktrace
#./bootstrap.sh --with-libraries=chrono,timer,program_options,filesystem,system,stacktrace
#./b2 --prefix=?? install ???



# SRC DIRS
make_dirs=deps/jellyfish-2.2.0 deps/seqan-library-2.0.0 deps/boost lib src tests

if MAKE_DOCS
    make_dirs += doc
endif

SUBDIRS = $(make_dirs)

AUTOMAKE_OPTIONS = subdir-objects