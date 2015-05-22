##########################################################################
# Setup
##########################################################################

# Autoconf initialistion.  Sets package name version and contact details
AC_INIT([Kmer Analysis Toolkit (KAT)],[1.1.0],[bernardo.clavijo@tgac.ac.uk and/or daniel.mapleson@tgac.ac.uk],[kat],[http://www.tgac.ac.uk/kat])

# Require autoconf 2.53 or higher
AC_PREREQ([2.53])

# Check that the configure script is being run from the right directory
AC_CONFIG_SRCDIR([src/kat.cc])

# Put auxillary build scripts in 'build-aux' directory
AC_CONFIG_AUX_DIR([build-aux])

# Set the macro directory
AC_CONFIG_MACRO_DIR([m4])

# Automake initialisation.  Sets some default flags.  We use 'foreign' because we have a markdown README file.
# I've only tried this with automake 1.11 and above, so you require at least that version installed on your 
# machine if you want to regenerate 'the Makefile.in's.  If you're looking at this comment in the 'configure' 
# script then this the am__api_version, will be set to whatever I used to create this file.  Or it might be 
# whatever you used if you made it yourself with autoconf (or autoreconf, or the autogen.sh script)
AM_INIT_AUTOMAKE([1.11 foreign -Wall no-define])

# Most of KAT is written in C++.  This is what this macro means: "Do compilation tests using CXX and CXXCPP and use extension .C for test programs. Use compilation flags: CPPFLAGS with CXXCPP, and both CPPFLAGS and CXXFLAGS with CXX."
AC_LANG_CPLUSPLUS

# Set output variable CXXCPP to a command that runs the C++ preprocessor. If ‘$CXX -E’ doesn't work, /lib/cpp is used. It is portable to run CXXCPP only on files with a .c, .C, .cc, or .cpp extension.
AC_PROG_CXXCPP

# AC_PROG_CXX determine a C++ compiler to use. Check whether the environment variable CXX or CCC (in that order) is set; if so, then set output variable CXX to its value.  We ignore the CXXFLAGS already set in the environment. This is necessary to avoid AC_PROG_XX automatically adding "-g -O2" to CXXFLAGS
CXXFLAGS=""
AC_PROG_CXX

# Determine a C compiler to use. If CC is not already set in the environment, check for gcc and cc, then for other C compilers. Set output variable CC to the name of the compiler found.
# NOTE: This probably isn't necessary.  I don't think we have any C code in KAT.
CFLAGS=""
AC_PROG_CC

# Ensure there's C++11 support
AX_CXX_COMPILE_STDCXX_11(,[mandatory])

# Create the 'install-sh' script, which is used as the 'install' make target
AC_PROG_INSTALL



######################################################################
# Debug compilation support.
######################################################################

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])

# If in debug mode make a really slow binary with debug symbols, otherwise make a fast one without debug symbols.  Seeing as these are being treated the same we could put them on CPPFLAGS instead, but I like to keep CPPFLAGS reserved for include paths and use CXX and C FLAGS for compiler settings.
if test x"$debugit" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="-g -Wno-uninitialized -O0"
    AM_CFLAGS="-g -Wno-uninitialized -O0"
    
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="-g -O3"
    AM_CFLAGS="-g -O3"
fi


######################################################################
# Dependency checks
######################################################################


# Check some (almost) standard functionality is present that we require to run KAT.
#AC_CHECK_HEADER_STDBOOL  # Commented this out because this was introduced after 2.63.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRERROR_R
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_SIZE_T
AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([sqrt])
AC_C_INLINE
AC_CHECK_HEADERS([string.h])
AC_CHECK_HEADERS([iostream])
AC_CHECK_HEADERS([stdio.h])
AC_CHECK_HEADERS([fstream])
AC_CHECK_HEADERS([vector])
AC_CHECK_HEADERS([getopt.h])



# Check pthreads exists
AX_PTHREAD(, [AC_MSG_ERROR([pthreads lib not found.  Please ensure that pthreads is properly built and configured.])])


## Check for boost
AX_BOOST_BASE([1.52],, [AC_MSG_ERROR([Boost not found.  Please ensure that boost is properly built and the BOOST_ROOT environment variable is set.  Alternatively you can override BOOST_ROOT with the --with-boost option.])])
AX_BOOST_FILESYSTEM
AX_BOOST_SYSTEM
AX_BOOST_TIMER
AX_BOOST_THREAD
AX_BOOST_UNIT_TEST_FRAMEWORK

# Combine BOOST variables (apart for BOOST_TEST)
BOOST_LDFLAGS="${BOOST_FILESYSTEM_LDFLAGS} ${BOOST_PROGRAM_OPTIONS_LDFLAGS} ${BOOST_SYSTEM_LDFLAGS} ${BOOST_TIMER_LDFLAGS} ${BOOST_THREAD_LDFLAGS}"
BOOST_LIBS="${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${BOOST_SYSTEM_LIB} ${BOOST_TIMER_LIB} ${BOOST_THREAD_LIB}"

AM_CXXFLAGS="-DCPLUSPLUS"
AC_SUBST([AM_CXXFLAGS])

AM_CPPFLAGS="${BOOST_CPPFLAGS}"
AC_SUBST([AM_CPPFLAGS])

AM_LDFLAGS="${BOOST_LDFLAGS}"
AC_SUBST([AM_LDFLAGS])

AM_LIBS="${BOOST_LIBS} ${PTHREAD_CFLAGS}"
AC_SUBST([AM_LIBS])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile deps/jellyfish-2.2.0/Makefile src/Makefile tests/Makefile])
AC_OUTPUT
