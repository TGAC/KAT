AM_LDFLAGS = -all-static -pthread
AM_CPPFLAGS = -Wall -Wno-sign-compare -Wnon-virtual-dtor -I../src/inc
AM_CXXFLAGS = -g

TESTS = check_sect check_comp
check_PROGRAMS = check_sect check_comp

check_sect_SOURCES = check_sect.cc ../src/sect/sect.hpp

check_comp_SOURCES = check_comp.cc ../src/comp/comp.hpp

check_sect_LDADD = ../lib/libz.a
