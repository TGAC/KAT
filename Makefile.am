
# Tells aclocal to put 'aclocal.m4' into the 'm4' directory
ACLOCAL_AMFLAGS = -I m4

# ADDITIONAL FILES TO INSTALL
EXTRA_DIST = include \
README.md \
COPYING \
AUTHORS \
deps/jellyfish-2.2.0/**.cc \
deps/jellyfish-2.2.0/**.cpp \
deps/jellyfish-2.2.0/**.hpp \
deps/jellyfish-2.2.0/**.h \
deps/jellyfish-2.2.0/**.i \
deps/jellyfish-2.2.0/**.py \
deps/jellyfish-2.2.0/**.pm \
deps/jellyfish-2.2.0/**.t \
deps/jellyfish-2.2.0/**.am \
deps/jellyfish-2.2.0/**.sh \
deps/jellyfish-2.2.0/**.in \
deps/jellyfish-2.2.0/**.m4 \
deps/jellyfish-2.2.0/doc/* \
deps/jellyfish-2.2.0/README \
deps/jellyfish-2.2.0/LICENSE \
deps/jellyfish-2.2.0/configure \
deps/jellyfish-2.2.0/config.guess \
deps/jellyfish-2.2.0/configure.ac \
deps/jellyfish-2.2.0/depcomp \
deps/jellyfish-2.2.0/gtest.mk \
deps/jellyfish-2.2.0/install-sh \
deps/jellyfish-2.2.0/missing \
deps/jellyfish-2.2.0/test-driver

# FILES NOT TO INSTALL
dist_noinst_SCRIPTS = autogen.sh clean.sh

# SCRIPTS TO ADD TO BIN DIR
dist_bin_SCRIPTS = \
        scripts/kat_comp_reads.sh \
        scripts/dist_analysis.py

# SRC DIRS
SUBDIRS = deps/jellyfish-2.2.0 src tests scripts


