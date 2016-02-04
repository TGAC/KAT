##########################################################################
# Setup
##########################################################################

# Autoconf initialistion.  Sets package name version and contact details
AC_PREREQ([2.68])
AC_INIT([kat],[2.1.0],[https://github.com/TGAC/KAT/issues],[kat],[https://github.com/TGAC/KAT])
AC_CONFIG_SRCDIR([src/kat.cc])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG_CPLUSPLUS
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_INSTALL

# Ensure there's C++11 support
AX_CXX_COMPILE_STDCXX_11(,[mandatory])

# Automake configuration.
AM_INIT_AUTOMAKE([1.11 subdir-objects foreign no-define tar-ustar])
AM_SILENT_RULES([yes])
m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR

# Libtool setup
LT_INIT

# Detect OS
UNAME_S=`uname -s`
if test x"${UNAME_S}" == x"Linux"; then
    AC_DEFINE([OS_LINUX], [1], [Linux detected])
fi
if test x"${UNAME_S}" == x"Darwin"; then
    AC_DEFINE([OS_MAC], [1], [Mac detected])
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
AC_CHECK_HEADERS([memory])
AC_CHECK_HEADERS([vector])
AC_CHECK_HEADERS([glob.h])


# Check pthreads exists
AX_PTHREAD(, [AC_MSG_ERROR([pthreads lib not found.  Please ensure that pthreads is properly built and configured.])])


# Plotting
AC_ARG_ENABLE([gnuplot], AS_HELP_STRING([--enable-gnuplot], [Enable gnuplot plotting even if python matplotlib is available]))

AM_PATH_PYTHON([], [], AC_MSG_ERROR([Python interpreter not found.]))

AX_PYTHON_DEVEL([>= '3.1'])

pymod_good="no"
if test x"${PYTHON_VERSION}" != x""; then
    if test x"${PYTHON_EXTRA_LIBS}" = x""; then
        pymod_good="no"
        AC_MSG_WARN([Python3 detected but Python3 development library was not found.  If you wish to use python plotting please install python3 library.  e.g. "sudo apt-get install python3-dev" on debian systems.])
    fi
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
fi

# Check for gnuplot on path.  Just emit a warning if not present
if test x"${pymod_good}" = x"no" || test "x$enable_gnuplot" = "xyes"; then
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

# Check for sphinx to build documentation (optional)
AC_CHECK_PROG([sphinx], [sphinx-build], [yes], [no])
if test x"${sphinx}" = x"yes"; then
	#Now check if the version is correct
	sphinx_major_version=`sphinx-build --version 2>&1 | cut -d' ' -f3 | cut -d '.' -f 1 | tr -d " \t\n\r"`
        sphinx_minor_version=`sphinx-build --version 2>&1 | cut -d' ' -f3 | cut -d '.' -f 2 | tr -d " \t\n\r"`
	if test x"${sphinx_major_version}" = x""; then sphinx_major_version=0; fi
	if test x"${sphinx_minor_version}" = x""; then sphinx_minor_version=0; fi
	major_size=${#sphinx_major_version}
	minor_size=${#sphinx_minor_version}
	if (( ("${major_size}" > 2) || ("${minor_size}" > 2) )); then sphinx_major_version=0; sphinx_minor_version=0; fi
	if (( ("${sphinx_major_version}" < 1 ) || (  ( "${sphinx_major_version}" >= 1 ) && ( "${sphinx_minor_version}" < 3 ) ) )); then
    		AC_MSG_WARN([Sphinx detected, but version is too old.  Wont create documentation.  You can still find the documentation online at: https://kat.readthedocs.org/en/latest/])
		sphinx="no"
	fi
else
    AC_MSG_WARN([Sphinx not detected, wont create documentation.  You can still find the documentation online at: https://kat.readthedocs.org/en/latest/])    
fi
AM_CONDITIONAL([MAKE_DOCS], [test x$sphinx = xyes])
AC_SUBST([MAKE_DOCS])


## Check for boost
AX_BOOST_BASE([1.52],, [AC_MSG_ERROR([Boost not found.  Please ensure that boost is properly built and the BOOST_ROOT environment variable is set.  Alternatively you can override BOOST_ROOT with the --with-boost option.])])
AX_BOOST_FILESYSTEM
AX_BOOST_SYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_CHRONO
AX_BOOST_TIMER


define([PC_FILE], lib/kat-2.1.pc)


# Combine BOOST variables (apart for BOOST_TEST)
BOOST_LIBS="${BOOST_FILESYSTEM_LIB} ${BOOST_PROGRAM_OPTIONS_LIB} ${BOOST_SYSTEM_LIB} ${BOOST_TIMER_LIB}"
BOOST_STATIC_LIBS="${BOOST_TIMER_STATIC_LIB} ${BOOST_CHRONO_STATIC_LIB} ${BOOST_FILESYSTEM_STATIC_LIB} ${BOOST_PROGRAM_OPTIONS_STATIC_LIB} ${BOOST_SYSTEM_STATIC_LIB} "
AC_SUBST([BOOST_STATIC_LIBS])

AM_CXXFLAGS="-DCPLUSPLUS"
AC_SUBST([AM_CXXFLAGS])

if test x"${pymod_good}" = x"yes"; then 
    AM_CPPFLAGS="${BOOST_CPPFLAGS} ${PYTHON_CPPFLAGS}"
    AM_LIBS="${PTHREAD_CFLAGS} ${BOOST_STATIC_LIBS} ${PYTHON_EXTRA_LIBS} ${PYTHON_LIBS}"
    AM_LDFLAGS="${BOOST_LDFLAGS} ${PYTHON_EXTRA_LDFLAGS}"
else
    AM_CPPFLAGS="${BOOST_CPPFLAGS}"
    AM_LIBS="${PTHREAD_CFLAGS} ${BOOST_STATIC_LIBS}"
    AM_LDFLAGS="${BOOST_LDFLAGS}"
fi

AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([AM_LIBS])


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile doc/Makefile doc/source/conf.py lib/kat-2.1.pc lib/Makefile src/Makefile tests/Makefile tests/compat.sh deps/seqan-library-2.0.0/Makefile])
AC_CONFIG_SUBDIRS([deps/jellyfish-2.2.0])
AC_OUTPUT

if test x"${pymod_good}" = x"yes" && test "x$enable_gnuplot" = "xno"; then
    AC_MSG_NOTICE([Using python plotting])
else
    if test x"${gnuplot_found}" = x"yes"; then
        if test "x$enable_gnuplot" = "xno"; then
            AC_MSG_WARN([Python3, or some required python modules where not found.  Python via matplotlib is the preferred plotting method in KAT.  If you want to use python plotting (the preferred option) then we suggest you install anaconda3 (or python3 with numpy, scipy and matplotlib separately).  However, gnuplot was detected so KAT plots will be generated via gnuplot instead.])
        else
            AC_MSG_NOTICE([Using gnuplot plotting])
        fi
    else
        AC_MSG_WARN([Neither python (with plotting modules present) or Gnuplot not detected.  KAT will still work (minus plots). Should you require plotting functionality, you must install either anaconda3 (or python3 with numpy, scipy and matplotlib) or gnuplot, and then reconfigure and compile to enable plotting in KAT.])
    fi
fi 

if test x"${sphinx}" = x"no"; then
    AC_MSG_WARN([Sphinx not detected, or version is too old.  Wont create documentation.  You can still find the documentation online at: https://kat.readthedocs.org/en/latest/])    
fi


