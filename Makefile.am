## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS  = -I m4

# Installed artifacts
dist_doc_DATA             = README.rst
nobase_pkginclude_HEADERS = tuna/tuna.h
lib_LTLIBRARIES           = tuna/libtuna.la

# The sources build into the library
tuna_libtuna_la_SOURCES = tuna/tuna.c

# Example programs
noinst_PROGRAMS         =
noinst_PROGRAMS        += examples/basic
examples_basic_SOURCES  = basic/basic.c
examples_basic_LDADD    = tuna/libtuna.la

# TODO Tests

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
