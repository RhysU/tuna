## Process this file with automake to produce Makefile.in
## Notice non-recursive Automake employed
ACLOCAL_AMFLAGS  = -I m4

# Installed artifacts
dist_doc_DATA             = README.rst
nobase_pkginclude_HEADERS = tuna/tuna.h
nobase_pkginclude_HEADERS = tuna/stats.h
lib_LTLIBRARIES           = tuna/libtuna.la

# The sources build into the library
tuna_libtuna_la_SOURCES  = tuna/tuna.c  tuna/tuna.h
tuna_libtuna_la_SOURCES += tuna/stats.c tuna/stats.h

# Example programs
noinst_PROGRAMS         =
noinst_PROGRAMS        += examples/basic
examples_basic_SOURCES  = examples/basic.c
examples_basic_LDADD    = tuna/libtuna.la

# TODO Tests

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
