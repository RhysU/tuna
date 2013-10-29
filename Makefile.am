## Process this file with automake to produce Makefile.in
## Notice non-recursive Automake employed
ACLOCAL_AMFLAGS  = -I m4
noinst_PROGRAMS  =

# Installed artifacts
dist_doc_DATA              = README.rst
nobase_pkginclude_HEADERS  = tuna/ltqnorm.h
nobase_pkginclude_HEADERS += tuna/stats.h
nobase_pkginclude_HEADERS += tuna/tuna.h
lib_LTLIBRARIES            = tuna/libtuna.la

# The sources build into the library
tuna_libtuna_la_SOURCES  = tuna/stats.c   tuna/stats.h
tuna_libtuna_la_SOURCES += tuna/tuna.c    tuna/tuna.h
tuna_libtuna_la_SOURCES += tuna/ltqnorm.c tuna/ltqnorm.h
tuna_libtuna_la_LIBADD   = $(LIBM)

# A minimum viable proof of concept
noinst_PROGRAMS      += examples/poc
examples_poc_SOURCES  = examples/poc.c
examples_poc_LDADD    = tuna/libtuna.la $(LIBM)

# Example programs
noinst_PROGRAMS        += examples/basic
examples_basic_SOURCES  = examples/basic.c
examples_basic_LDADD    = tuna/libtuna.la $(LIBM)

# Test programs
TESTS                = tests/stats
noinst_PROGRAMS     += tests/stats
tests_stats_SOURCES  = tests/stats.c tests/fct.h
tests_stats_LDADD    = tuna/libtuna.la $(LIBM)

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
