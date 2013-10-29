AC_INIT([tuna], [0.0.1], [rhys.ulerich@gmail.com])
AC_CONFIG_SRCDIR([tuna/tuna.h])
AC_CONFIG_HEADERS([tuna/config.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([1.11 check-news color-tests foreign parallel-tests silent-rules subdir-objects tar-ustar -Wall])
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C99
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_SEARCH_LIBS([clock_gettime],[rt])
LT_INIT
LT_LIB_M
AC_SUBST([LIBTOOL_DEPS])
AC_CONFIG_FILES([
  Makefile
])
AC_OUTPUT
