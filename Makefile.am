##############################################################################
# Automake configuration
##############################################################################

AM_CXXFLAGS =  -I. -Wall -Wextra -pedantic -Wno-c++11-extensions -std=c++11 
AM_CXXFLAGS += -g -O3 $(GCC_SPECIFIC_CXXFLAGS)
AM_LDFLAGS  =  -no-undefined

ACLOCAL_AMFLAGS	= -I m4

## libsemigroups headers
pkginclude_HEADERS =  src/libsemigroups-debug.h 
pkginclude_HEADERS += src/blocks.h               src/cong.h 
pkginclude_HEADERS += src/elements.h             src/semigroups.h 
pkginclude_HEADERS += src/rws.h                  src/rwse.h 
pkginclude_HEADERS += src/semiring.h             src/partition.h
pkginclude_HEADERS += src/recvec.h               src/report.h	
pkginclude_HEADERS += src/timer.h                src/uf.h

nodist_include_HEADERS = config/libsemigroups-config.h

lib_LTLIBRARIES = libsemigroups.la

## libsemigroups sources

libsemigroups_la_SOURCES =  src/blocks.cc   src/cong.cc
libsemigroups_la_SOURCES += src/elements.cc src/semigroups.cc
libsemigroups_la_SOURCES += src/rws.cc	    src/rwse.cc
libsemigroups_la_SOURCES += src/uf.cc

libsemigroups_la_CPPFLAGS = -DCONFIG_H $(CODE_COVERAGE_CPPFLAGS)
libsemigroups_la_CFLAGS = $(CODE_COVERAGE_CFLAGS) 
libsemigroups_la_CXXFLAGS = $(AM_CXXFLAGS) $(CODE_COVERAGE_CXXFLAGS)

libsemigroups_la_LIBS = $(AM_LDFLAGS) $(CODE_COVERAGE_LIBS)

## Extra files for the distribution

LINT_FORMAT_FILES_EXTRA =  src/cong/kbfp.h src/cong/kbp.h 
LINT_FORMAT_FILES_EXTRA += src/cong/p.h    src/cong/tc.h

LINT_FORMAT_FILES_EXTRA += src/cong/kbfp.cc src/cong/kbp.cc 
LINT_FORMAT_FILES_EXTRA += src/cong/p.cc    src/cong/tc.cc

EXTRA_DIST = $(LINT_FORMAT_FILES_EXTRA)
EXTRA_DIST += tests/catch.hpp
EXTRA_DIST += README.md LICENSE CPPLINT.cfg .clang-format Doxyfile

BENCHMARK_LINT_FORMAT =  benchmark/src/cong.bm.cpp
BENCHMARK_LINT_FORMAT += benchmark/src/semigroups.bm.cpp
BENCHMARK_LINT_FORMAT += benchmark/src/nridempotents.bm.cpp
BENCHMARK_LINT_FORMAT += benchmark/src/examples.h

## lstest sources 

CHECK_LOG_DIR = tests/log
CHECK_PROG    = ./lstest$(EXEEXT)
CHECK_FLAGS1  = -d yes --force-colour --order lex --abort 
CHECK_FLAGS2  = [quick],[standard]

check_PROGRAMS = lstest

lstest_SOURCES =  tests/main.test.cc
lstest_SOURCES += tests/blocks.test.cc	  tests/recvec.test.cc
lstest_SOURCES += tests/cong.test.cc	  tests/rws.test.cc
lstest_SOURCES += tests/elements.test.cc  tests/rwse.test.cc
lstest_SOURCES += tests/kbfp.test.cc	  tests/semigroups.test.cc
lstest_SOURCES += tests/kbp.test.cc       tests/semiring.test.cc
lstest_SOURCES += tests/p.test.cc         tests/tc.test.cc
lstest_SOURCES += tests/partition.test.cc tests/uf.test.cc

lstest_CPPFLAGS = -DCONFIG_H $(CODE_COVERAGE_CPPFLAGS)
lstest_CFLAGS = $(CODE_COVERAGE_CFLAGS) 
lstest_CXXFLAGS = $(AM_CXXFLAGS) $(CODE_COVERAGE_CXXFLAGS)

lstest_LIBS = $(AM_LDFLAGS) $(CODE_COVERAGE_LIBS)
lstest_LDADD = libsemigroups.la

## Make targets

@CODE_COVERAGE_RULES@

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	rm -rf $(CHECK_LOG_DIR) $(check_PROGRAMS)
	rm -rf benchmark/src/build
	rm -f benchmark/gbench/*.pyc
	rm -f benchmark/gbench/__py_cache__/*.pyc

superclean: distclean
	git clean -xdf --exclude *.swp --exclude *.swo
.PHONY: superclean

check:	
	@$(mkdir_p) $(CHECK_LOG_DIR)
	@today=`date +%Y-%m-%d-%H-%M-%S`;                                   \
	echo "Running the tests ($(CHECK_LOG_DIR)/$$today.log) . . .";      \
	$(CHECK_PROG) $(CHECK_FLAGS1) $(CHECK_FLAGS2) | tee -a $(CHECK_LOG_DIR)/$$today.log; \
	( ! grep -q -E "FAILED|failed" $(CHECK_LOG_DIR)/$$today.log )

check-quick: CHECK_FLAGS2 = [quick]
check-quick: check

check-all: CHECK_FLAGS2 = 
check-all: check
.PHONY: check-quick, check-all

format: 
	@echo "Running clang-format . . .";  
	@clang-format -i $(pkginclude_HEADERS) $(libsemigroups_la_SOURCES) 
	@clang-format -i $(LINT_FORMAT_FILES_EXTRA) $(lstest_SOURCES) 
	@clang-format -i $(BENCHMARK_LINT_FORMAT)

lint: 
	@echo "Running cpplint . . .";  
	@cpplint $(pkginclude_HEADERS) $(libsemigroups_la_SOURCES) 
	@cpplint $(LINT_FORMAT_FILES_EXTRA) $(lstest_SOURCES) 	
	@cpplint $(BENCHMARK_LINT_FORMAT)
	
doc:
	@echo "Generating documentation . . ."; \
	doxygen Doxyfile
	@echo "See: html/index.html";
