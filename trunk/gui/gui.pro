TEMPLATE = app
QT = gui \
 core \
 xml \
 svg \
 dbus
CONFIG += qt \
 warn_on \
 console \
 release \
 qdbus
LIBS = -lgslcblas -lgsl
SOURCES = mainwindowimpl.cpp \
 gui.cpp \
 plotimpl.cpp \
 gelmanrubinimpl.cpp \
 paramqt.cpp \
 ../warg/src/rng.cpp \
 ../warg/src/wargxml.cpp \
 ../warg/src/param.cpp \
 ../warg/src/rectree.cpp \
 ../warg/src/tree.cpp \
 ../warg/src/node.cpp \
 ../warg/src/data.cpp \
 ../warg/src/recedge.cpp \
 paramcons.cpp \
 paramconsmult.cpp \
 densityontree.cpp \
 outputfile.cpp \
 parammr.cpp \
 paramtreecons.cpp \
 ../warg/src/metropolis.cpp \
 ../warg/src/move.cpp \
 ../warg/src/movetheta.cpp \
 ../warg/src/moveremedge.cpp \
 ../warg/src/moveaddedge.cpp \
 ../warg/src/movesitechange.cpp \
 ../warg/src/movetimechange.cpp \
 ../warg/src/moverho.cpp \
 ../warg/src/moveedgechange.cpp \
 ../warg/src/movedelta.cpp \
 ../warg/src/moveageclonal.cpp \
 ../warg/src/movescaletree.cpp \
 ../warg/src/moveregraftclonal.cpp \
 ../warg/src/movegreedytree.cpp \
 ../warg/src/mpiutils.cpp \
 colouredimpl.cpp \
 pdimpl.cpp \
 pheatimpl.cpp \
 heatimpl.cpp
HEADERS = mainwindowimpl.h \
 mainapplicationdbusadaptor.h \
 plotimpl.h \
 gelmanrubinimpl.h \
 paramqt.h \
 ../warg/src/rng.h \
 ../warg/src/wargxml.h \
 ../warg/src/param.h \
 ../warg/src/rectree.h \
 ../warg/src/tree.h \
 ../warg/src/node.h \
 ../warg/src/rng.h \
 ../warg/src/data.h \
 ../warg/src/recedge.h \
 paramcons.h \
 paramconsmult.h \
 densityontree.h \
 outputfile.h \
 parammr.h \
 paramtreecons.h \
 ../warg/src/metropolis.h \
 ../warg/src/move.h \
 ../warg/src/movetheta.h \
 ../warg/src/moveremedge.h \
 ../warg/src/moveaddedge.h \
 ../warg/src/movesitechange.h \
 ../warg/src/movetimechange.h \
 ../warg/src/moverho.h \
 ../warg/src/moveedgechange.h \
 ../warg/src/movedelta.h \
 ../warg/src/moveageclonal.h \
 ../warg/src/movescaletree.h \
 ../warg/src/moveregraftclonal.h \
 ../warg/src/movegreedytree.h \
 ../warg/src/mpiutils.h \
 colouredimpl.h \
 pdimpl.h \
 pheatimpl.h \
 heatimpl.h
FORMS += mainwindow.ui \
 plot.ui \
 gelmanrubin.ui \
 coloured.ui \
 pd.ui \
 pheat.ui \
 heat.ui
