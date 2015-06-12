TEMPLATE = app
CONFIG += release
LIBS = -Wl,--no-as-needed -lgslcblas -lgsl
QMAKE_CXXFLAGS = -g -funroll-loops
SOURCES = src/data.cpp \
 src/rng.cpp \
 src/wargxml.cpp \
 src/move.cpp \
 src/metropolis.cpp \
 src/node.cpp \
 src/param.cpp \
 src/recedge.cpp \
 src/tree.cpp \
 src/rectree.cpp \
 src/weakarg.cpp \
 src/movetheta.cpp \
 src/moveremedge.cpp \
 src/moveaddedge.cpp \
 src/movesitechange.cpp \
 src/movetimechange.cpp \
 src/moverho.cpp \
 src/moveedgechange.cpp \
 src/movedelta.cpp \
 src/moveageclonal.cpp \
 src/movescaletree.cpp \
 src/moveregraftclonal.cpp \
 src/movegreedytree.cpp \
 src/mpiutils.cpp
HEADERS = src/data.h \
 src/rng.h \
 src/wargxml.h \
 src/move.h \
 src/metropolis.h \
 src/node.h \
 src/param.h \
 src/recedge.h \
 src/tree.h \
 src/rectree.h \
 src/movetheta.h \
 src/rng.h \
 src/moveremedge.h \
 src/moveaddedge.h \
 src/movesitechange.h \
 src/movetimechange.h \
 src/moverho.h \
 src/moveedgechange.h \
 src/movedelta.h \
 src/moveageclonal.h \
 src/movescaletree.h \
 src/moveregraftclonal.h \
 src/movegreedytree.h \
 src/mpiutils.h
