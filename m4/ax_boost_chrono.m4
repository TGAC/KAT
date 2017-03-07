# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_boost_chrono.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_CHRONO
#
# DESCRIPTION
#
#   Test for System library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_CHRONO_LIB)
#     AC_SUBST(BOOST_CHRONO_STATIC_LIB)
#
#   And sets:
#
#     HAVE_BOOST_CHRONO
#
# LICENSE
#
#   Copyright (c) 2012 Xiyue Deng <manphiz@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_BOOST_CHRONO],
[
    AC_ARG_WITH([boost-chrono],
    AS_HELP_STRING([--with-boost-chrono@<:@=special-lib@:>@],
               [use the Chrono library from boost - it is possible to specify a certain library for the linker
                    e.g. --with-boost-chrono=boost_chrono-gcc-mt ]),
    [
    if test "$withval" = "no"; then
        want_boost="no"
    elif test "$withval" = "yes"; then
        want_boost="yes"
        ax_boost_user_chrono_lib=""
    else
        want_boost="yes"
        ax_boost_user_chrono_lib="$withval"
        fi
    ],
    [want_boost="yes"]
    )

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
        AC_REQUIRE([AC_CANONICAL_BUILD])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::Chrono library is available,
					   ax_cv_boost_chrono,
            [AC_LANG_PUSH([C++])
                CXXFLAGS_SAVE=$CXXFLAGS

                AC_COMPILE_IFELSE(  [AC_LANG_PROGRAM([[@%:@include <boost/chrono.hpp>]],
                                [[boost::chrono::system_clock::time_point time;]])],
                                ax_cv_boost_chrono=yes, ax_cv_boost_chrono=no)
                                CXXFLAGS=$CXXFLAGS_SAVE
                                AC_LANG_POP([C++])
            ]
        )
        
        if test "x$ax_cv_boost_chrono" = "xyes"; then
            AC_SUBST(BOOST_CPPFLAGS)

            AC_DEFINE(HAVE_BOOST_CHRONO,,[define if the Boost::Chrono library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`

            ax_lib=""
            ax_static_lib=""
            no_find="yes"
            no_link="yes"
            link_chrono="no"
            link_chrono_static="no"

            LDFLAGS_SAVE=$LDFLAGS
            if test "x$ax_boost_user_chrono_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_chrono*.so* $BOOSTLIBDIR/libboost_chrono*.dylib* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_chrono.*\)\.so.*$;\1;' -e 's;^lib\(boost_chrono.*\)\.dylib.*$;\1;'` ; do
                    ax_lib=${libextension}
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_CHRONO_LIB="-l$ax_lib"; AC_SUBST(BOOST_CHRONO_LIB) link_chrono="yes"; break],
                        [link_chrono="no"])
                done
                if test "x$link_chrono" != "xyes"; then
                    for libextension in `ls $BOOSTLIBDIR/boost_chrono*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_chrono.*\)\.dll.*$;\1;'` ; do
                        ax_lib=${libextension}
                        AC_CHECK_LIB($ax_lib, exit,
                            [BOOST_CHRONO_LIB="-l$ax_lib"; AC_SUBST(BOOST_CHRONO_LIB) link_chrono="yes"; break],
                            [link_chrono="no"])
                    done
                fi
                for libextension in `ls $BOOSTLIBDIR/libboost_chrono*.a* 2>/dev/null` ; do
                    ax_static_lib=${libextension}
                    AC_CHECK_FILE($ax_static_lib,
                        [BOOST_CHRONO_STATIC_LIB="$ax_static_lib"; AC_SUBST(BOOST_CHRONO_STATIC_LIB) link_chrono_static="yes"; break],
                        [link_chrono_static="no"])
                done

                no_find="no"
                if [[ -z "$ax_lib" ]] && [[ -z "$ax_static_lib" ]]; then
                    no_find="yes"
                fi

                no_link="no"
                if [[ "$link_chrono" == "no" ]] && [[ "$link_chrono_static" == "no" ]]; then
                    no_link="yes"
                fi

            else
                for ax_lib in $ax_boost_user_chrono_lib boost_chrono-$ax_boost_user_chrono_lib; do
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_CHRONO_LIB="-l$ax_lib"; AC_SUBST(BOOST_CHRONO_LIB) link_chrono="yes"; break],
                        [link_chrono="no"])
                done

            fi
            if [[ -z "$ax_lib" ]]; then
                AC_MSG_WARN(Could not find a dynamic version of boost_chrono)
            fi
            if [[ -z "$ax_static_lib" ]]; then
                AC_MSG_WARN(Could not find a static version of boost_chrono)
            fi
            if [[ "$no_find" == "yes" ]]; then
                AC_MSG_ERROR(Could not find any version boost_chrono to link to)
            fi

            if [[ "$link_chrono" == "no" ]]; then
                AC_MSG_WARN(Could not dynamic link against boost_chrono)
            fi
            if [[ "$link_chrono_static" == "no" ]]; then
                AC_MSG_WARN(Could not static link against boost_chrono)
            fi
            if [[ "$no_link" == "yes" ]]; then
                AC_MSG_ERROR(Could not link against any boost_chrono lib)
            fi

        fi

        CPPFLAGS="$CPPFLAGS_SAVED"
        LDFLAGS="$LDFLAGS_SAVED"
    fi
])