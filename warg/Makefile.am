ACLOCAL_AMFLAGS = -I m4

.PHONY: profile
profile:
	make clean
	make CXXFLAGS="$(CXXFLAGS) -fprofile-generate"
	make check
	make mostlyclean
	make CXXFLAGS="$(CXXFLAGS) -fprofile-use"

CLEANFILES = src/*.gcda src/*.gcno


SUBDIRS = src scripts
EXTRA_DIST = weakarg.pro
