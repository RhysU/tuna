## Process this file with automake to produce Makefile.in
## Notice non-recursive Automake employed
ACLOCAL_AMFLAGS    = -I m4
AM_CPPFLAGS        = -I$(top_srcdir)/tuna -I$(top_srcdir)/cdflib
check_PROGRAMS     =
dist_check_SCRIPTS =
noinst_PROGRAMS    =
TESTS              =

# Installed artifacts
dist_doc_DATA    = README.rst
include_HEADERS  = tuna/tuna.h
lib_LTLIBRARIES  = tuna/libtuna.la

# Two small, C89-compliant source files comprise the entire library
# They can be copied directly into other projects if desired
tuna_libtuna_la_SOURCES = tuna/tuna.c tuna/tuna.h
tuna_libtuna_la_CFLAGS  = $(AM_CFLAGS) $(STDC89)
tuna_libtuna_la_LIBADD  = $(LIBM)

# C99-based example programs
noinst_PROGRAMS        += examples/basic
examples_basic_SOURCES  = examples/basic.c
examples_basic_CFLAGS   = $(AM_CFLAGS) $(STDC99)
examples_basic_LDADD    = tuna/libtuna.la $(LIBM)

noinst_PROGRAMS           += examples/algodive
examples_algodive_SOURCES  = examples/algodive.c
examples_algodive_CFLAGS   = $(AM_CFLAGS) $(STDC99)
examples_algodive_LDADD    = tuna/libtuna.la $(LIBM)

noinst_PROGRAMS            += examples/smallsort
examples_smallsort_SOURCES  = examples/smallsort.c
examples_smallsort_CFLAGS   = $(AM_CFLAGS) $(STDC99)
examples_smallsort_LDADD    = tuna/libtuna.la $(LIBM)

# C99-based test programs
# (either dedicated sources or simply reusing examples as tests)
TESTS               += tests/stats
noinst_PROGRAMS     += tests/stats
tests_stats_SOURCES  = tests/stats.c tests/fct.h
tests_stats_CFLAGS   = $(AM_CFLAGS) $(STDC99)
tests_stats_LDADD    = tuna/libtuna.la $(LIBM)

TESTS               += tests/welch
noinst_PROGRAMS     += tests/welch
tests_welch_SOURCES  = tests/welch.c tests/fct.h
tests_welch_CFLAGS   = $(AM_CFLAGS) $(STDC99)
tests_welch_LDADD    = tuna/libtuna.la $(LIBM)

TESTS               += tests/chunk
noinst_PROGRAMS     += tests/chunk
tests_chunk_SOURCES  = tests/chunk.c tests/fct.h
tests_chunk_CFLAGS   = $(AM_CFLAGS) $(STDC99)
tests_chunk_LDADD    = tuna/libtuna.la $(LIBM)

check_PROGRAMS       += examples/basic
dist_check_SCRIPTS   += tests/basic.sh
TESTS                += tests/basic.sh

check_PROGRAMS       += examples/smallsort
dist_check_SCRIPTS   += tests/smallsort.sh
TESTS                += tests/smallsort.sh

# CDFLIB might be handy for development but it is not installed
# include_HEADERS           = cdflib/cdflib.h
noinst_LTLIBRARIES          = cdflib/libcdflib.la
cdflib_libcdflib_la_SOURCES = cdflib/cdflib.c  cdflib/cdflib.h

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
