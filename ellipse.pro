TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

CONFIG += c++11


SOURCES += \
        algorithm/ellipse.c \
        algorithm/equation.c \
        algorithm/matrix.c \
        fileread.cpp \
        main.cpp

HEADERS += \
    algoheader.h \
    algorithm/algoheader.h \
    algorithm/ellipse.h \
    algorithm/equation.h \
    algorithm/matrix.h \
    fileread.h

