AC_INIT([tuna], [0.2.2], [rhys.ulerich@gmail.com])
AC_CONFIG_SRCDIR([tuna/tuna.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.11 check-news color-tests foreign parallel-tests silent-rules subdir-objects tar-ustar -Wall])
AC_PROG_CC
AM_PROG_CC_C_O
AX_CFLAGS_WARN_ALL
dnl Isolate -std=c89 flags into STDC89
CCSAVED=${CC}
AC_PROG_CC_C89
AX_CFLAGS_FORCE_C89([CC])
AC_SUBST([STDC89],[${CC#$CCSAVED}])
CC=$CCSAVED
dnl Isolate -std=c99 flags into STDC99
CCSAVED=${CC}
AC_PROG_CC_C99
AC_SUBST([STDC99],[${CC#$CCSAVED}])
CC=$CCSAVED
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_SEARCH_LIBS([clock_gettime],[rt])
AC_OPENMP
LT_INIT
LT_LIB_M
AC_SUBST([LIBTOOL_DEPS])
DX_DOXYGEN_FEATURE([ON])
DX_DOT_FEATURE([ON])
DX_HTML_FEATURE([ON])
DX_MAN_FEATURE([ON])
dnl Get version from git for Doxygen documentation
AC_MSG_CHECKING([for git describe version])
GIT_VERSION=`git describe --always --dirty 2>/dev/null`
if test -n "$GIT_VERSION"; then
    AC_MSG_RESULT([$GIT_VERSION])
else
    AC_MSG_RESULT([not available, using $PACKAGE_VERSION])
    GIT_VERSION=$PACKAGE_VERSION
fi
AC_SUBST([GIT_VERSION])
DX_INIT_DOXYGEN([tuna],[doxygen.cfg],[docs])
dnl Override VERSION with GIT_VERSION for Doxygen
DX_ENV_APPEND([VERSION], [$GIT_VERSION])
AC_CONFIG_FILES([
  Makefile
])
AC_OUTPUT
