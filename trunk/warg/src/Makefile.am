bin_PROGRAMS = warg

#OPTIMIZATION = -O1 -g -Wall
OPTIMIZATION = -O3 -Wall -mfpmath=sse -msse -msse2 -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations
warg_CXXFLAGS = $(GSL_CFLAGS) $(OPTIMIZATION) -Wall
warg_LDFLAGS = $(OPTIMIZATION)
warg_LDADD = $(GSL_LIBS)
warg_SOURCES = data.h rng.h wargxml.h exponential.h move.h metropolis.h moveaddedge.h moveageclonal.h movedelta.h moveedgechange.h moveregraftclonal.h moveremedge.h moverho.h movescaletree.h  movegreedytree.h movesitechange.h movetheta.h movetimechange.h mpiutils.h node.h param.h recedge.h rectree.h rng.h slotallocator.h tree.h weakarg.h data.cpp rng.cpp wargxml.cpp move.cpp metropolis.cpp moveaddedge.cpp moveageclonal.cpp movedelta.cpp moveedgechange.cpp moveregraftclonal.cpp moveremedge.cpp moverho.cpp movescaletree.cpp movegreedytree.cpp movesitechange.cpp movetheta.cpp movetimechange.cpp mpiutils.cpp node.cpp param.cpp recedge.cpp rectree.cpp tree.cpp weakarg.cpp 

check_SCRIPTS = profilerun.sh
TESTS = $(check_SCRIPTS)


profilerun.sh:
	echo './warg -s 532772982 -x 100 -y 50 testout' > profilerun.sh
	chmod +x profilerun.sh
