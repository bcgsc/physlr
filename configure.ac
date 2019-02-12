#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(PHYLSR, 0.0.1, cjustin@bcgsc.ca)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([overlap/MinimizerOverlap.cpp])
AC_CONFIG_HEADER([config.h])

#see http://stackoverflow.com/questions/8215785/how-can-i-use-ac-revision-with-git
AC_REVISION([m4_esyscmd([./tools/configure.commit])])
AC_DEFINE([GIT_REVISION], "[m4_esyscmd([./tools/configure.commit])]", [git revision information])

AC_PROG_RANLIB

# Checks for programs.
# make sure -O3 by default for speed
if test -z $CXXFLAGS; then
    AC_SUBST(CPPFLAGS,"-O3 $CPPFLAGS")
fi

AC_SUBST(CPPFLAGS,"-std=c++11 $CPPFLAGS")
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Werror $AM_CXXFLAGS")

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
#AC_CHECK_LIB([], [])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Options to configure.
AC_ARG_WITH(sparsehash, AS_HELP_STRING([--with-sparsehash=PATH],
			[specify directory for the google sparsehash header files]))
if test "$with_sparsehash" -a -d "$with_sparsehash"; then
	sparsehash_cppflags="-isystem$with_sparsehash"
fi

# Find the absolute paths.
srcdir=$(cd $srcdir; pwd)

# Set compiler flags.
AC_SUBST(CPPFLAGS,"-I$srcdir $sparsehash_cppflags $CPPFLAGS")
#AC_SUBST(LDFLAGS, "$LDFLAGS")

# Checks for header files.
AC_LANG([C++])	#sets the auto conf to look for c++ headers
#AC_CHECK_HEADERS([zlib.h])

AC_CHECK_HEADERS([google/dense_hash_map])
if test $ac_cv_header_google_dense_hash_map != yes; then
	AC_MSG_ERROR([Requires google sparsehash])
fi

# Check for OpenMP.
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	OPENMP_CXXFLAGS=-Wno-unknown-pragmas
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

# Checks for typedefs, structures, and compiler characteristics.

AC_CONFIG_FILES([
		Makefile
		overlap/Makefile
])
AC_OUTPUT

