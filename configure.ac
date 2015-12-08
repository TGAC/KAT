##########################################################################
# Setup
##########################################################################

# Autoconf initialistion.  Sets package name version and contact details
AC_INIT([kat],[2.1.0],[bernardo.clavijo@tgac.ac.uk and/or daniel.mapleson@tgac.ac.uk],[kat],[http://www.tgac.ac.uk/kat])

# Require autoconf 2.53 or higher
AC_PREREQ([2.53])

# Check that the configure script is being run from the right directory
AC_CONFIG_SRCDIR([src/kat.cc])

# Put auxillary build scripts in 'build-aux' directory
AC_CONFIG_AUX_DIR([build-aux])

# Set the macro directory
AC_CONFIG_MACRO_DIR([m4])

# Automake initialisation.
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign no-define tar-ustar])

# Use this to display simplified compilation output
#AM_SILENT_RULES([yes])

# Most of KAT is written in C++.  This is what this macro means: "Do compilation tests using CXX and CXXCPP and use extension .C for test programs. Use compilation flags: CPPFLAGS with CXXCPP, and both CPPFLAGS and CXXFLAGS with CXX."
AC_LANG_CPLUSPLUS

# Set output variable CXXCPP to a command that runs the C++ preprocessor. If ‘$CXX -E’ doesn't work, /lib/cpp is used. It is portable to run CXXCPP only on files with a .c, .C, .cc, or .cpp extension.
AC_PROG_CXXCPP

# Identify a c++ compiler to use
AC_PROG_CXX

# Ensure there's C++11 support
AX_CXX_COMPILE_STDCXX_11(,[mandatory])

# Create the 'install-sh' script, which is used as the 'install' make target
AC_PROG_INSTALL

# Required for using the archiver
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR


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
AC_CHECK_HEADERS([memory])
AC_CHECK_HEADERS([vector])
AC_CHECK_HEADERS([glob.h])


# Check pthreads exists
AX_PTHREAD(, [AC_MSG_ERROR([pthreads lib not found.  Please ensure that pthreads is properly built and configured.])])

AX_PYTHON_DEVEL([>='3.2.0'])

pymod_good="yes"
AX_PYTHON_MODULE(numpy, ,python3)
if test x"${PYMOD}" = x"no"; then 
	pymod_good="no" 
fi
AX_PYTHON_MODULE(matplotlib, ,python3)
if test x"${PYMOD}" = x"no"; then 
	pymod_good="no" 
fi
AX_PYTHON_MODULE(scipy, ,python3)
if test x"${PYMOD}" = x"no"; then 
	pymod_good="no" 
fi

# Check for gnuplot on path.  Just emit a warning if not present
if test x"${pymod_good}" = x"no"; then
    AC_DEFINE([HAVE_PYTHON], [0], [Python not present])
    AC_CHECK_PROG([gnuplot_found], [gnuplot], [yes], [no])
    if test x"${gnuplot_found}" = x"no"; then
        AC_DEFINE([HAVE_GNUPLOT], [0], [Gnuplot not present])
    else
        AC_DEFINE([HAVE_GNUPLOT], [1], [Gnuplot present])
    fi
else
    AC_DEFINE([HAVE_PYTHON], [1], [Python present])
fi


## Check for boost
AX_BOOST_BASE([1.52],, [AC_MSG_ERROR([Boost not found.  Please ensure that boost is properly built and the BOOST_ROOT environment variable is set.  Alternatively you can override BOOST_ROOT with the --with-boost option.])])
AX_BOOST_FILESYSTEM
AX_BOOST_SYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_CHRONO
AX_BOOST_TIMER
AX_BOOST_UNIT_TEST_FRAMEWORK


# Combine BOOST variables (apart for BOOST_TEST)
BOOST_LIBS="${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${BOOST_SYSTEM_LIB} ${BOOST_TIMER_LIB}"
BOOST_STATIC_LIBS="${BOOST_TIMER_STATIC_LIB} ${BOOST_CHRONO_STATIC_LIB} ${BOOST_FILESYSTEM_STATIC_LIB} ${BOOST_PROGRAM_OPTIONS_STATIC_LIB} ${BOOST_SYSTEM_STATIC_LIB} "
AC_SUBST([BOOST_STATIC_LIBS])

AM_CXXFLAGS="-DCPLUSPLUS"
AC_SUBST([AM_CXXFLAGS])

AM_LDFLAGS="${BOOST_LDFLAGS} -Xlinker -export-dynamic"
AC_SUBST([AM_LDFLAGS])

if test x"${pymod_good}" = x"yes"; then 
    AM_CPPFLAGS="${BOOST_CPPFLAGS} -I${PYTHON_INCLUDE_DIR}"
    AM_LIBS="${PTHREAD_CFLAGS} ${BOOST_STATIC_LIBS} -ldl -lutil -lm -l${PYTHON_LIB}"
else
    AM_CPPFLAGS="${BOOST_CPPFLAGS}"
    AM_LIBS="${PTHREAD_CFLAGS} ${BOOST_STATIC_LIBS}"
fi

AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LIBS])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile deps/seqan-library-2.0.0/Makefile])
AC_CONFIG_SUBDIRS([deps/jellyfish-2.2.0])
AC_OUTPUT

if test x"${pymod_good}" = x"yes"; then
    echo "Using python plotting"
else
    if test x"${gnuplot_found}" = x"yes"; then
        AC_MSG_WARN([Python3, or some required python modules where not found.  Python via matplotlib is the preferred plotting method in KAT.  If you want to use python plotting (the preferred option) then we suggest you install anaconda3 (or python3 with numpy, scipy and matplotlib separately).  However, gnuplot was detected so KAT plots will be generated via gnuplot instead.])
    else
        AC_MSG_WARN([Neither python (with plotting modules present) or Gnuplot not detected.  KAT will still work (minus plots). Should you require plotting functionality, you must install either anaconda3 (or python3 with numpy, scipy and matplotlib) or gnuplot, and then reconfigure and compile to enable plotting in KAT.])
    fi
fi 


