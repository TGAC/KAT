AC_INIT([Kmer Analysis Toolkit (KAT)], [0.3.5], [daniel.mapleson@tgac.ac.uk], [kat], [http://www.tgac.ac.uk])
AM_INIT_AUTOMAKE([1.11 -Wall no-define])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG_CPLUSPLUS

# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
# check for C++ preprocessor and compiler and the library compiler
# (might change the compiler flags)
AC_PROG_CXXCPP
AC_PROG_CXX
# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS



AC_PROG_INSTALL

##########################################################################
# debug compilation support
##########################################################################

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])

if test x"$debugit" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wno-uninitialized -O0"
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -O3"
fi

##########################################################################


######################################################################
# DOXYGEN SUPPORT
######################################################################

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN([$PACKAGE_NAME],[doxygen.cfg])

######################################################################


######################################################################
# Configure External Libraries
######################################################################


# Allow users to specify custom directories for jellyfish, boost and seqan

seqan_dir=":"
boost_dir=":"
doxygen_dir=":"


PKG_CHECK_MODULES([JELLYFISH], [jellyfish-1.1 >= 1.1.10])

AC_ARG_WITH(seqan,
        [  --with-seqan=LIB_DIR    seqan root directory is here. ],
        [ seqan_dir=$withval], [ seqan_dir=":"])

AC_ARG_WITH(boost,
        [  --with-boost=LIB_DIR    boost root directory is here... (boost is only required for running unit tests) ],
        [ boost_dir=$withval], [ boost_dir=":"])

AC_ARG_WITH(doxygen,
        [  --with-doxygen=LIB_DIR  doxygen root directory is here... (doxygen is only required for generating source code documentation) ],
        [ doxygen_dir=$withval], [ doxygen_dir=":"])




# Modify environment variables as appropriate

if test "${boost_dir}" != ":"; then
    CPPFLAGS="$CPPFLAGS -I${boost_dir}/include"
    LIBS="$LIBS ${boost_dir}/lib/libboost_unit_test_framework.so"
fi

if test "${seqan_dir}" != ":"; then
    CPPFLAGS="$CPPFLAGS -I${seqan_dir}/include"
fi

if test "${doxygen_dir}" != ":"; then
    PATH="$PATH ${doxygen_dir}/bin"
fi


# Only jellyfish and seqan are required. boost is only used for unit testing
# if the tests fail it's not the end of the world, and doxygen is only required
# for generating source code documentation.


AC_CHECK_HEADER([seqan/sequence.h], [], [
         echo "Seqan is required to compile and run KAT.  Please make sure seqan is installed on your system."
         echo "If seqan is installed but the headers are not currently available through the CPPFLAGS environment"
         echo "variable then please re-run \"./configure\" with the \"--with-seqan\" argument to specify the root"
         echo "directory of your seqan installation"
         exit -1],
         [])



# distribute the changed variables among the Makefiles

AC_SUBST([JELLYFISH_CFLAGS])
AC_SUBST([JELLYFISH_LIBS])
AC_SUBST([LIBS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
