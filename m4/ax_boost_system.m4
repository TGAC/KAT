# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_boost_system.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_SYSTEM
#
# DESCRIPTION
#
#   Test for System library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_SYSTEM_LIB)
#     AC_SUBST(BOOST_SYSTEM_STATIC_LIB)
#
#   And sets:
#
#     HAVE_BOOST_SYSTEM
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Michael Tindal
#   Copyright (c) 2008 Daniel Casimiro <dan.casimiro@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 17

AC_DEFUN([AX_BOOST_SYSTEM],
[
	AC_ARG_WITH([boost-system],
	AS_HELP_STRING([--with-boost-system@<:@=special-lib@:>@],
                   [use the System library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-system=boost_system-gcc-mt ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_system_lib=""
        else
		    want_boost="yes"
		ax_boost_user_system_lib="$withval"
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

        AC_CACHE_CHECK(whether the Boost::System library is available,
					   ax_cv_boost_system,
        [AC_LANG_PUSH([C++])
			 CXXFLAGS_SAVE=$CXXFLAGS

			 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/system/error_code.hpp>]],
                                   [[boost::system::system_category]])],
                   ax_cv_boost_system=yes, ax_cv_boost_system=no)
			 CXXFLAGS=$CXXFLAGS_SAVE
             AC_LANG_POP([C++])
		])
		        if test "x$ax_cv_boost_system" = "xyes"; then
            AC_SUBST(BOOST_CPPFLAGS)

            AC_DEFINE(HAVE_BOOST_SYSTEM,,[define if the Boost::System library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`

            ax_lib=""
            ax_static_lib=""
            no_find="yes"
            no_link="yes"
            link_system="no"
            link_system_static="no"

            LDFLAGS_SAVE=$LDFLAGS
            if test "x$ax_boost_user_system_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_system*.so* $BOOSTLIBDIR/libboost_system*.dylib* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_system.*\)\.so.*$;\1;' -e 's;^lib\(boost_system.*\)\.dylib.*$;\1;'` ; do
                    ax_lib=${libextension}
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_SYSTEM_LIB="-l$ax_lib"; AC_SUBST(BOOST_SYSTEM_LIB) link_system="yes"; break],
                        [link_system="no"])
                done
                if test "x$link_system" != "xyes"; then
                    for libextension in `ls $BOOSTLIBDIR/boost_system*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_system.*\)\.dll.*$;\1;'` ; do
                        ax_lib=${libextension}
                        AC_CHECK_LIB($ax_lib, exit,
                            [BOOST_SYSTEM_LIB="-l$ax_lib"; AC_SUBST(BOOST_SYSTEM_LIB) link_system="yes"; break],
                            [link_system="no"])
                    done
                fi
                for libextension in `ls $BOOSTLIBDIR/libboost_system*.a* 2>/dev/null` ; do
                    ax_static_lib=${libextension}
                    AC_CHECK_FILE($ax_static_lib,
                        [BOOST_SYSTEM_STATIC_LIB="$ax_static_lib"; AC_SUBST(BOOST_SYSTEM_STATIC_LIB) link_system_static="yes"; break],
                        [link_system_static="no"])
                done

                no_find="no"
                if [[ -z "$ax_lib" ]] && [[ -z "$ax_static_lib" ]]; then
                    no_find="yes"
                fi

                no_link="no"
                if [[ "$link_system" == "no" ]] && [[ "$link_system_static" == "no" ]]; then
                    no_link="yes"
                fi

            else
                for ax_lib in $ax_boost_user_system_lib boost_system-$ax_boost_user_system_lib; do
                    AC_CHECK_LIB($ax_lib, exit,
                        [BOOST_SYSTEM_LIB="-l$ax_lib"; AC_SUBST(BOOST_SYSTEM_LIB) link_system="yes"; break],
                        [link_system="no"])
                done

            fi
            if [[ -z "$ax_lib" ]]; then
                AC_MSG_WARN(Could not find a dynamic version of boost_system)
            fi
            if [[ -z "$ax_static_lib" ]]; then
                AC_MSG_WARN(Could not find a static version of boost_system)
            fi
            if [[ "$no_find" == "yes" ]]; then
                AC_MSG_ERROR(Could not find any version boost_system to link to)
            fi

            if [[ "$link_system" == "no" ]]; then
                AC_MSG_WARN(Could not dynamic link against boost_system)
            fi
            if [[ "$link_system_static" == "no" ]]; then
                AC_MSG_WARN(Could not static link against boost_system)
            fi
            if [[ "$no_link" == "yes" ]]; then
                AC_MSG_ERROR(Could not link against any boost_system lib)
            fi

        fi

	CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
    fi
])
