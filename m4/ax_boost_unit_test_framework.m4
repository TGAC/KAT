# ================================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_boost_unit_test_framework.html
# ================================================================================
#
# SYNOPSIS
#
#   AX_BOOST_UNIT_TEST_FRAMEWORK
#
# DESCRIPTION
#
#   Test for Unit_Test_Framework library from the Boost C++ libraries. The
#   macro requires a preceding call to AX_BOOST_BASE. Further documentation
#   is available at <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIB)
#     AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_STATIC_LIB)
#
#   And sets:
#
#     HAVE_BOOST_UNIT_TEST_FRAMEWORK
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 19

AC_DEFUN([AX_BOOST_UNIT_TEST_FRAMEWORK],
[
	AC_ARG_WITH([boost-unit-test-framework],
	AS_HELP_STRING([--with-boost-unit-test-framework@<:@=special-lib@:>@],
                   [use the Unit_Test_Framework library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-unit-test-framework=boost_unit_test_framework-gcc ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_unit_test_framework_lib=""
        else
		    want_boost="yes"
		ax_boost_user_unit_test_framework_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::Unit_Test_Framework library is available,
					   ax_cv_boost_unit_test_framework,
            [AC_LANG_PUSH([C++])
                             AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/test/unit_test.hpp>]],
                                        [[using boost::unit_test::test_suite;
                                                             test_suite* test= BOOST_TEST_SUITE( "Unit test example 1" ); return 0;]])],
                       ax_cv_boost_unit_test_framework=yes, ax_cv_boost_unit_test_framework=no)
             AC_LANG_POP([C++])
            ]
        )
		
        if test "x$ax_cv_boost_unit_test_framework" = "xyes"; then
            AC_DEFINE(HAVE_BOOST_UNIT_TEST_FRAMEWORK,,[define if the Boost::Unit_Test_Framework library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
            LDFLAGS_SAVE=$LDFLAGS
            if test "x$ax_boost_user_unit_test_framework_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_unit_test_framework*.so* $BOOSTLIBDIR/libboost_unit_test_framework*.dylib* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_unit_test_framework.*\)\.so.*$;\1;' -e 's;^lib\(boost_unit_test_framework.*\)\.dylib.*$;\1;'` ; do
                    ax_lib=${libextension}
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_UNIT_TEST_FRAMEWORK_LIB="-l$ax_lib"; AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIB) link_unit_test_framework="yes"; break],
                        [link_unit_test_framework="no"])
                done
                if test "x$link_unit_test_framework" != "xyes"; then
                    for libextension in `ls $BOOSTLIBDIR/boost_unit_test_framework*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_unit_test_framework.*\)\.dll.*$;\1;'` ; do
                        ax_lib=${libextension}
                        AC_CHECK_LIB($ax_lib, exit,
                            [BOOST_UNIT_TEST_FRAMEWORK_LIB="-l$ax_lib"; AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIB) link_unit_test_framework="yes"; break],
                            [link_unit_test_framework="no"])
                    done
                fi
                for libextension in `ls $BOOSTLIBDIR/libboost_unit_test_framework*.a* 2>/dev/null` ; do
                    ax_static_lib=${libextension}
                    AC_CHECK_FILE($ax_static_lib,
                        [BOOST_UNIT_TEST_FRAMEWORK_STATIC_LIB="$ax_static_lib"; AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_STATIC_LIB) link_unit_test_framework_static="yes"; break],
                        [link_unit_test_framework_static="no"])
                done

                no_find="no"
                if test "x$ax_lib" = "x"; then
                    if test "x$ax_static_lib" = "x"; then
                        no_find="yes"
                    fi
                fi

                no_link="no"
                if test "x$link_unit_test_framework" != "xyes"; then
                    if test "x$link_unit_test_framework_static" != "xyes"; then
                        no_link="yes"
                    fi
                fi

            else
                for ax_lib in $ax_boost_user_unit_test_framework_lib boost_unit_test_framework-$ax_boost_user_unit_test_framework_lib; do
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_UNIT_TEST_FRAMEWORK_LIB="-l$ax_lib"; AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIB) link_unit_test_framework="yes"; break],
                        [link_unit_test_framework="no"])
                done

            fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_WARN(Could not find a dynamic version of the library!)
            elif test "x$ax_static_lib" = "x"; then
                AC_MSG_WARN(Could not find a static version of the library!)
            fi
            if test "x$no_find" = "xyes"; then
                AC_MSG_ERROR(Could not find any version of the library to link to)
            fi

            if test "x$link_unit_test_framework" = "xno"; then
                AC_MSG_WARN(Could not dynamic link against $ax_lib !)
            elif test "x$link_unit_test_framework_static" = "xno"; then
                AC_MSG_WARN(Could not static link against $ax_static_lib!)
            fi
            if test "x$no_link" = "xyes"; then
                AC_MSG_ERROR(Could not link against any boost-unit_test_framework lib)
            fi

        fi

		CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
	fi
])