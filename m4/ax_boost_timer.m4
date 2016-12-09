# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_boost_timer.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_TIMER
#
# DESCRIPTION
#
#   Test for System library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_TIMER_LIB)
#     AC_SUBST(BOOST_TIMER_STATIC_LIB)
#
#   And sets:
#
#     HAVE_BOOST_TIMER
#
# LICENSE
#
#   Copyright (c) 2012 Xiyue Deng <manphiz@gmail.com>
#   Copyright (c) 2012 Murray Cumming <murrayc@openismus.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 2 (based on serial 1 of ax_boost_locale.m4 with some simple find/replace by Murray Cumming)

AC_DEFUN([AX_BOOST_TIMER],
[
	AC_ARG_WITH([boost-timer],
	AS_HELP_STRING([--with-boost-timer@<:@=special-lib@:>@],
                   [use the Timer library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-timer=boost_timer-gcc-mt ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_timer_lib=""
        else
		    want_boost="yes"
		ax_boost_user_timer_lib="$withval"
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

        AC_CACHE_CHECK(whether the Boost::Timer library is available,
					   ax_cv_boost_timer,
            [AC_LANG_PUSH([C++])
			 CXXFLAGS_SAVE=$CXXFLAGS

			 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/timer/timer.hpp>]],
                                   [[boost::timer::cpu_timer().stop();]])],
                   ax_cv_boost_timer=yes, ax_cv_boost_timer=no)
			 CXXFLAGS=$CXXFLAGS_SAVE
             AC_LANG_POP([C++])
            ]
        )
		
        if test "x$ax_cv_boost_timer" = "xyes"; then
            AC_SUBST(BOOST_CPPFLAGS)

            AC_DEFINE(HAVE_BOOST_TIMER,,[define if the Boost::Timer library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`

            LDFLAGS_SAVE=$LDFLAGS
            if test "x$ax_boost_user_timer_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_timer*.so* $BOOSTLIBDIR/libboost_timer*.dylib* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_timer.*\)\.so.*$;\1;' -e 's;^lib\(boost_timer.*\)\.dylib.*$;\1;'` ; do
                    ax_lib=${libextension}
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_TIMER_LIB="-l$ax_lib"; AC_SUBST(BOOST_TIMER_LIB) link_timer="yes"; break],
                        [link_timer="no"])
                done
                if test "x$link_timer" != "xyes"; then
                    for libextension in `ls $BOOSTLIBDIR/boost_timer*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_timer.*\)\.dll.*$;\1;'` ; do
                        ax_lib=${libextension}
                        AC_CHECK_LIB($ax_lib, exit,
                            [BOOST_TIMER_LIB="-l$ax_lib"; AC_SUBST(BOOST_TIMER_LIB) link_timer="yes"; break],
                            [link_timer="no"])
                    done
                fi
                for libextension in `ls $BOOSTLIBDIR/libboost_timer*.a* 2>/dev/null` ; do
                    ax_static_lib=${libextension}
                    AC_CHECK_FILE($ax_static_lib,
                        [BOOST_TIMER_STATIC_LIB="$ax_static_lib"; AC_SUBST(BOOST_TIMER_STATIC_LIB) link_timer_static="yes"; break],
                        [link_timer_static="no"])
                done

                no_find="no"
                if test "x$ax_lib" = "x"; then
                    if test "x$ax_static_lib" = "x"; then
                        no_find="yes"
                    fi
                fi

                no_link="no"
                if test "x$link_timer" != "xyes"; then
                    if test "x$link_timer_static" != "xyes"; then
                        no_link="yes"
                    fi
                fi

            else
                for ax_lib in $ax_boost_user_timer_lib boost_timer-$ax_boost_user_timer_lib; do
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_TIMER_LIB="-l$ax_lib"; AC_SUBST(BOOST_TIMER_LIB) link_timer="yes"; break],
                        [link_timer="no"])
                done

            fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_WARN(Could not find a dynamic version of boost_timer)
            elif test "x$ax_static_lib" = "x"; then
                AC_MSG_WARN(Could not find a static version of boost_timer)
            fi
            if test "x$no_find" = "xyes"; then
                AC_MSG_ERROR(Could not find any version of the library to link to)
            fi

            if test "x$link_timer" = "xno"; then
                AC_MSG_WARN(Could not dynamic link against $ax_lib)
            elif test "x$link_timer_static" = "xno"; then
                AC_MSG_WARN(Could not static link against $ax_static_lib)
            fi
            if test "x$no_link" = "xyes"; then
                AC_MSG_ERROR(Could not link against any boost_timer lib)
            fi

        fi
        
        CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
    fi
])