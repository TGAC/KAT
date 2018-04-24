
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
dist_noinst_SCRIPTS = autogen.sh antigen.sh build_boost.sh

# SRC DIRS
make_dirs=deps/jellyfish-2.2.0 deps/seqan-library-2.0.0 lib src tests

if MAKE_PYKAT
    make_dirs += scripts
endif

if MAKE_DOCS
    make_dirs += doc
endif

SUBDIRS = $(make_dirs)

AUTOMAKE_OPTIONS = subdir-objects
