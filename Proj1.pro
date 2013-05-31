TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    trialfct.cpp \
    Proj1.cpp \
    hamilton.cpp \
    mcint.cpp \
    positions.cpp \
    function.cpp \
    lib.cpp \
    trialfct_molecule.cpp

LIBS += -larmadillo -lconfig++

HEADERS += \
    trialfct.h \
    hamilton.h \
    mcint.h \
    positions.h \
    function.h \
    lib.h \
    trialfct_molecule.h
