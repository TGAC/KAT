
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST =	README.md \
		COPYING \
		AUTHORS \
		INSTALL \
		ChangeLog


# Scripts NOT to install
dist_noinst_SCRIPTS = autogen.sh clean.sh

# Scripts to install
dist_bin_SCRIPTS = scripts/dist_analysis.py

# SRC DIRS
SUBDIRS = deps/jellyfish-2.2.0 deps/seqan-library-2.0.0 src tests

AUTOMAKE_OPTIONS = subdir-objects
