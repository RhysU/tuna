## Process this file with automake to produce Makefile.in
## Notice non-recursive Automake employed
ACLOCAL_AMFLAGS  = -I m4
noinst_PROGRAMS  =
TESTS            =

# Installed artifacts
dist_doc_DATA              = README.rst
nobase_pkginclude_HEADERS  = tuna/algo.h
nobase_pkginclude_HEADERS += tuna/as3.h
nobase_pkginclude_HEADERS += tuna/cdflib.h
nobase_pkginclude_HEADERS += tuna/countof.h
nobase_pkginclude_HEADERS += tuna/kernel.h
nobase_pkginclude_HEADERS += tuna/ltqnorm.h
nobase_pkginclude_HEADERS += tuna/rand.h
nobase_pkginclude_HEADERS += tuna/stats.h
nobase_pkginclude_HEADERS += tuna/tuna.h
nobase_pkginclude_HEADERS += tuna/welch.h
lib_LTLIBRARIES            = tuna/libtuna.la

# The sources build into the library
tuna_libtuna_la_SOURCES  = tuna/algo.c    tuna/algo.h
tuna_libtuna_la_SOURCES += tuna/as3.c     tuna/as3.h
tuna_libtuna_la_SOURCES += tuna/cdflib.c  tuna/cdflib.h
tuna_libtuna_la_SOURCES += tuna/kernel.c  tuna/kernel.h
tuna_libtuna_la_SOURCES += tuna/ltqnorm.c tuna/ltqnorm.h
tuna_libtuna_la_SOURCES += tuna/rand.c    tuna/rand.h
tuna_libtuna_la_SOURCES += tuna/stats.c   tuna/stats.h
tuna_libtuna_la_SOURCES += tuna/tuna.c    tuna/tuna.h
tuna_libtuna_la_SOURCES += tuna/welch.c   tuna/welch.h
tuna_libtuna_la_LIBADD   = $(LIBM)

# Example programs
noinst_PROGRAMS      += examples/poc
examples_poc_SOURCES  = examples/poc.c
examples_poc_LDADD    = tuna/libtuna.la $(LIBM)

noinst_PROGRAMS        += examples/basic
examples_basic_SOURCES  = examples/basic.c
examples_basic_LDADD    = tuna/libtuna.la $(LIBM)

noinst_PROGRAMS           += examples/algodive
examples_algodive_SOURCES  = examples/algodive.c
examples_algodive_LDADD    = tuna/libtuna.la $(LIBM)

noinst_PROGRAMS            += examples/smallsort
examples_smallsort_SOURCES  = examples/smallsort.c
examples_smallsort_LDADD    = tuna/libtuna.la $(LIBM)

# Test programs
TESTS               += tests/stats
noinst_PROGRAMS     += tests/stats
tests_stats_SOURCES  = tests/stats.c tests/fct.h
tests_stats_LDADD    = tuna/libtuna.la $(LIBM)

TESTS               += tests/welch
noinst_PROGRAMS     += tests/welch
tests_welch_SOURCES  = tests/welch.c tests/fct.h
tests_welch_LDADD    = tuna/libtuna.la $(LIBM)

TESTS                += tests/kernel
noinst_PROGRAMS      += tests/kernel
tests_kernel_SOURCES  = tests/kernel.c tests/fct.h
tests_kernel_LDADD    = tuna/libtuna.la $(LIBM)

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool
