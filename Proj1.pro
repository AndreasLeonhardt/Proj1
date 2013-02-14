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
    hamilton_analytical.cpp \
    hamilton_numerical.cpp

LIBS += -larmadillo -lconfig++

HEADERS += \
    trialfct.h \
    hamilton.h \
    mcint.h \
    positions.h \
    hamilton_analytical.h \
    hamilton_numerical.h
