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


# check for jellyfish boost and seqan, which are dependencies for kat
jellyfish_dir=":"
boost_dir=":"
seqan_dir=":"

AC_ARG_WITH(jellyfish,
        [  --with-jellyfish=LIB_DIR     search static jellyfish library there. ],
        [ jellyfish_dir=$withval], [ jellyfish_dir=":"])


AC_ARG_WITH(boost,
        [  --with-boost=LIB_DIR     search static boost library there. ],
        [ boost_dir=$withval], [ boost_dir=":"])


AC_ARG_WITH(seqan,
        [  --with-boost=LIB_DIR     search static seqan library there. ],
        [ seqan_dir=$withval], [ seqan_dir=":"])



if test "$jellyfish_dir" = ":"; then
    echo "WARNING: JELLYFISH: Did not specify jellyfish v1.1.10 directory... assuming env vars are set correctly"

else
    echo
    echo "JELLYFISH: Custom jellyfish dir specified.  Modifying environment variables for jellyfish: ${jellyfish_dir}"

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


if test "${boost_dir}" = ":"; then
    echo "WARNING: BOOST: Did not specify boost directory... assuming env vars are set correctly"

else
    echo
    echo "BOOST: Custom boost dir specified.  Modifying environment variables for boost: ${boost_dir}"

    CPPFLAGS="$CPPFLAGS -I${boost_dir}/include"
    LIBS="$LIBS ${boost_dir}/lib/libboost_unit_test_framework.so"

    echo
    echo "+++++++++++++++++++++++++++++++++++++"
    echo "CFLAGS now: ${CFLAGS}"
    echo "CPPFLAGS now: ${CPPFLAGS}"
    echo "CXXFLAGS now: ${CXXFLAGS}"
    echo "LDFLAGS now: ${LDFLAGS}"
    echo "LIBS now: ${LIBS}"
    echo "+++++++++++++++++++++++++++++++++++++"
fi


if test "${seqan_dir}" = ":"; then
    echo "WARNING: SEQAN: Did not specify seqan directory... assuming env vars are set correctly"

else
    echo
    echo "SEQAN: Custom seqan dir specified.  Modifying environment variables for seqan: ${seqan_dir}"

    CPPFLAGS="$CPPFLAGS -I${seqan_dir}/include"

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
