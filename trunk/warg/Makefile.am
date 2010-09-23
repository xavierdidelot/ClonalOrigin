ACLOCAL_AMFLAGS = -I m4

.PHONY: profile
profile:
	make clean
	make CXXFLAGS="$(CXXFLAGS) -fprofile-generate"
	make check
	make mostlyclean
	make CXXFLAGS="$(CXXFLAGS) -fprofile-use"

CLEANFILES = src/*.gcda src/*.gcno


SUBDIRS = src
EXTRA_DIST = convert.pl rescale.pl weakarg.pro
