## Process this file with automake to produce Makefile.in
## Notice non-recursive Automake employed
ACLOCAL_AMFLAGS    = -I m4
AM_CPPFLAGS        = -I$(top_srcdir)/tuna
check_PROGRAMS     =
dist_check_SCRIPTS =
dist_doc_DATA      =
EXTRA_DIST         =
include_HEADERS    =
lib_LTLIBRARIES    =
noinst_LTLIBRARIES =
noinst_PROGRAMS    =
TESTS              =

# Build a debug library for internal, performance agnostic use...
# (Two small, C89-compliant source files comprise the entire library)
# (These files can be copied directly into other projects if desired)
noinst_LTLIBRARIES             += tuna/libtuna-debug.la
tuna_libtuna_debug_la_SOURCES   = tuna/tuna.c tuna/tuna.h
tuna_libtuna_debug_la_CFLAGS    = $(AM_CFLAGS) $(STDC89)
tuna_libtuna_debug_la_CPPFLAGS  = $(AM_CPPFLAGS)
tuna_libtuna_debug_la_LIBADD    = $(LIBM)

# ...and reusing the libtuna-debug settings, prepare an assertless library.
dist_doc_DATA            += README.rst
include_HEADERS          += tuna/tuna.h
lib_LTLIBRARIES          += tuna/libtuna.la
tuna_libtuna_la_SOURCES   = $(tuna_libtuna_debug_la_SOURCES)
tuna_libtuna_la_CFLAGS    = $(tuna_libtuna_debug_la_CFLAGS)
tuna_libtuna_la_CPPFLAGS  = $(tuna_libtuna_debug_la_CPPFLAGS) -DNDEBUG
tuna_libtuna_la_LIBADD    = $(tuna_libtuna_debug_la_LIBADD)
tuna_libtuna_la_LDFLAGS   = -release $(PACKAGE_VERSION)

# C99-based example programs with sources installed alongside the documentation
# Use libtuna-debug.la when performance is noncritical to improve coverage.
noinst_PROGRAMS        += examples/basic
examples_basic_SOURCES  = examples/basic.c
examples_basic_CFLAGS   = $(AM_CFLAGS) $(STDC99)
examples_basic_LDADD    = tuna/libtuna-debug.la $(LIBM)
dist_doc_DATA          += $(examples_basic_SOURCES)

noinst_PROGRAMS           += examples/algodive
examples_algodive_SOURCES  = examples/algodive.c
examples_algodive_CFLAGS   = $(AM_CFLAGS) $(STDC99)
examples_algodive_LDADD    = tuna/libtuna-debug.la $(LIBM)
dist_doc_DATA             += $(examples_algodive_SOURCES)

noinst_PROGRAMS            += examples/smallsort
examples_smallsort_SOURCES  = examples/smallsort.c
examples_smallsort_CFLAGS   = $(AM_CFLAGS) $(STDC99)
examples_smallsort_LDADD    = tuna/libtuna.la $(LIBM)
dist_doc_DATA              += $(examples_smallsort_SOURCES)

noinst_PROGRAMS            += examples/blockedmm
examples_blockedmm_SOURCES  = examples/blockedmm.c
examples_blockedmm_CFLAGS   = $(AM_CFLAGS) $(STDC99) $(OPENMP_CFLAGS)
examples_blockedmm_LDADD    = tuna/libtuna.la $(LIBM) $(OPENMP_CFLAGS)
dist_doc_DATA              += $(examples_blockedmm_SOURCES)

# C99-based test programs
# (either dedicated sources or simply reusing examples as tests)
TESTS               += tests/stats
noinst_PROGRAMS     += tests/stats
tests_stats_SOURCES  = tests/stats.c tests/fct.h
tests_stats_CFLAGS   = $(AM_CFLAGS) $(STDC99)
tests_stats_LDADD    = tuna/libtuna-debug.la $(LIBM)

TESTS               += tests/welch
noinst_PROGRAMS     += tests/welch
tests_welch_SOURCES  = tests/welch.c tests/fct.h
tests_welch_CFLAGS   = $(AM_CFLAGS) $(STDC99)
tests_welch_LDADD    = tuna/libtuna-debug.la $(LIBM)

TESTS               += tests/chunk
noinst_PROGRAMS     += tests/chunk
tests_chunk_SOURCES  = tests/chunk.c tests/fct.h
tests_chunk_CFLAGS   = $(AM_CFLAGS) $(STDC99)
tests_chunk_LDADD    = tuna/libtuna-debug.la $(LIBM)

check_PROGRAMS       += examples/basic
dist_check_SCRIPTS   += tests/basic.sh
TESTS                += tests/basic.sh

check_PROGRAMS       += examples/smallsort
dist_check_SCRIPTS   += tests/smallsort.sh
TESTS                += tests/smallsort.sh

check_PROGRAMS       += examples/blockedmm
dist_check_SCRIPTS   += tests/blockedmm.sh
TESTS                += tests/blockedmm.sh

# Ease the process of running Artistic Style on the library sources
EXTRA_DIST += astyle.rc
astyle:     astyle.out
astyle.out: $(tuna_libtuna_la_SOURCES)
	astyle --verbose --options=$(srcdir)/astyle.rc $? | tee astyle.out

# Automates updating libtool dependencies when necessary
LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool

# Documentation hooks in with the 'html' target
html-local: doxygen-run

# Target for all documentation
docs: doxygen-run

# Required for AX_PROG_DOXYGEN
include doxygen.am
