##########################################################################
# Setup
##########################################################################

# Autoconf initialistion.  Sets package name version and contact details
AC_PREREQ([2.68])
AC_INIT([kat],[2.4.0],[https://github.com/TGAC/KAT/issues],[kat],[https://github.com/TGAC/KAT])
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
#Check for zlib
AX_CHECK_ZLIB(, [AC_MSG_ERROR([zlib not found.  Please ensure that zlib is properly built and configured.])])

# Required for Boost chrono on linux (clock_gettime).  Doesn't exist on Mac
if test x"${UNAME_S}" != x"Darwin"; then
    AC_SEARCH_LIBS([clock_gettime], [rt],
        [RT_LIB="-lrt"],
        [AC_MSG_ERROR([rt lib not found.  Please ensure that rt is properly built and configured.])])
else
    RT_LIB=""
fi

AC_ARG_ENABLE([pykat],
    AS_HELP_STRING([--disable-pykat], [Disable build of python package pykat.  If this is disabled, KAT will not create plots or perform extra distribution analysis.]), make_pykat="no", make_pykat="yes")

if [[ "${make_pykat}" == "yes" ]]; then

    AX_PYTHON_DEVEL([>= '3.5'])

    if [[ -z "${PYTHON_VERSION}" ]]; then
        AC_MSG_ERROR([Python3 not detected. If running on a debian system try "sudo apt-get install python3-dev".])
    fi

    if [[ -z "${PYTHON_EXTRA_LIBS}" ]]; then
        AC_MSG_ERROR([Python3 detected but Python3 development library was not found. If running on a debian system try "sudo apt-get install python3-dev".])
    fi

    pybin=python${PYTHON_VERSION}

    # Check for required python libraries
    AX_PYTHON_MODULE(numpy, [AC_MSG_ERROR([numpy not found])], ${pybin})
    AX_PYTHON_MODULE(matplotlib, [AC_MSG_ERROR([matplotlib not found])], ${pybin})
    AX_PYTHON_MODULE(scipy, [AC_MSG_ERROR([scipy not found])], ${pybin})

    #AC_DEFINE([HAVE_PYTHON], [1], [Python present])
    #AC_DEFINE_UNQUOTED([KAT_SITE_PKGS], ["${libdir}/$pybin/site-packages"], [KAT module site-packages at installed location])
fi

AM_CONDITIONAL([MAKE_PYKAT], [test x$make_pykat = xyes])
AC_SUBST([MAKE_PYKAT])

# Check for sphinx to build documentation (optional)
AC_CHECK_PROG([sphinx], [sphinx-build], [yes], [no])
if [[ "${sphinx}" == "yes" ]]; then
	#Now check if the version is correct
	sphinx_major_version=`sphinx-build --version 2>&1 | cut -d' ' -f3 | cut -d '.' -f 1 | tr -d " \t\n\r"`
        sphinx_minor_version=`sphinx-build --version 2>&1 | cut -d' ' -f3 | cut -d '.' -f 2 | tr -d " \t\n\r"`
	if [[ -z "${sphinx_major_version}" ]] ; then sphinx_major_version=0; fi
	if [[ -z "${sphinx_minor_version}" ]] ; then sphinx_minor_version=0; fi
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


define([PC_FILE], lib/kat-2.2.pc)

AM_CXXFLAGS=""
if test x"${UNAME_S}" == x"Linux"; then
	AM_CXXFLAGS="-static-libstdc++"
fi

AC_SUBST([AM_CXXFLAGS])

if [[ "${make_pykat}" == "yes" ]]; then
    AM_CPPFLAGS="${ZLIB_CPPFLAGS} ${PYTHON_CPPFLAGS}"
    AM_LIBS="${PTHREAD_CFLAGS} ${ZLIB_LIB} ${RT_LIB} ${PYTHON_EXTRA_LIBS} ${PYTHON_LIBS}"
    AM_LDFLAGS="${ZLIB_LDFLAGS} ${PYTHON_EXTRA_LDFLAGS} ${PYTHON_RPATH}"
else
    AM_CPPFLAGS="${ZLIB_CPPFLAGS}"
    AM_LIBS="${PTHREAD_CFLAGS} ${ZLIB_LIB} ${RT_LIB}"
    AM_LDFLAGS="${ZLIB_LDFLAGS}"
fi

AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([AM_LIBS])


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile doc/Makefile doc/source/conf.py lib/kat.pc lib/Makefile src/Makefile tests/Makefile tests/compat.sh deps/seqan-library-2.0.0/Makefile scripts/setup.py scripts/kat/__init__.py scripts/Makefile])
AC_CONFIG_SUBDIRS([deps/jellyfish-2.2.0])
AC_OUTPUT

if [[ "${make_pykat}" == "no" ]]; then
    AC_MSG_NOTICE([Python functionality disabled.  No plotting or distribution analysis will be performed with this KAT build.])
else
    AC_MSG_NOTICE([Python configuration is compatible with KAT.  Plotting and distribution analyses will be performed by this build.])
fi

if [[ "${sphinx}" == "no" ]]; then
    AC_MSG_WARN([Sphinx not detected, or version is too old.  Wont create documentation.  You can still find the documentation online at: https://kat.readthedocs.org/en/latest/])
else
    AC_MSG_NOTICE([Sphinx detected.  Documentation will be built.])
fi
