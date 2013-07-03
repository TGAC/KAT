AC_INIT([Kmer Analysis Toolkit (KAT)], [0.1], [daniel.mapleson@tgac.ac.uk], [kat], [http://www.tgac.ac.uk])
AM_INIT_AUTOMAKE([1.11 -Wall no-define])
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL


echo
echo "+++++++++++++++++++++++++++++++++++++"
echo "CFLAGS now: ${CFLAGS}"
echo "CPPFLAGS now: ${CPPFLAGS}"
echo "CXXFLAGS now: ${CXXFLAGS}"
echo "LDFLAGS now: ${LDFLAGS}"
echo "LIBS now: ${LIBS}"
echo "+++++++++++++++++++++++++++++++++++++"


# check for jellyfish and zlib
jellyfish_dir=":"

AC_ARG_WITH(jellyfish-dir,
        [  --with-jellyfish-dir=LIB_DIR     search static jellyfish library there. ],
        [ jellyfish_dir=$withval], [ jellyfish_dir=""])

if test "${jellyfish_dir}" = ":"; then
    echo "Did not specify jellyfish directory... assuming env vars are set correctly"

else
    echo
    echo "JELLYFISH: Modifying environment variables for jellyfish: ${jellyfish_dir}"

    CPPFLAGS="$CPPFLAGS -I${jellyfish_dir}/include/jellyfish-1.1.10"
    LIBS="$LIBS ${jellyfish_dir}/lib/libjellyfish-1.1.a"

    echo
    echo "+++++++++++++++++++++++++++++++++++++"
    echo "CFLAGS now: ${CFLAGS}"
    echo "CPPFLAGS now: ${CPPFLAGS}"
    echo "CXXFLAGS now: ${CXXFLAGS}"
    echo "LDFLAGS now: ${LDFLAGS}"
    echo "LIBS now: ${LIBS}"
    echo "+++++++++++++++++++++++++++++++++++++"
fi



AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
