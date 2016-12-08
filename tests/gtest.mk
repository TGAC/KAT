##############################
# Gtest build.
##############################
# Build rules for libraries.
check_LTLIBRARIES = libgtest.la libgtest_main.la

libgtest_la_SOURCES = gtest/src/gtest-all.cc
libgtest_main_la_SOURCES = gtest/src/gtest_main.cc
libgtest_main_la_LIBADD = libgtest.la
libgtest_la_CPPFLAGS = -I$(srcdir)
libgtest_main_la_CPPFLAGS = -I$(srcdir)

GTEST_SRC = gtest/src/gtest-all.cc	\
	    gtest/src/gtest_main.cc	\
	    gtest/gtest.h


