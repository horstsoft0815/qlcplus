include(../../../variables.pri)

TEMPLATE = app
LANGUAGE = C++
TARGET   = aboutbox_test

QT      += testlib gui script
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

INCLUDEPATH += ../../src
INCLUDEPATH += ../../../engine/src
DEPENDPATH  += ../../src

QMAKE_LIBDIR += ../../src
QMAKE_LIBDIR += ../../../engine/src
#QMAKE_LIBDIR += ../../../engine/audio/plugins/btrack
LIBS         += -lqlcplusui -lqlcplusengine

# Test sources
SOURCES += aboutbox_test.cpp
HEADERS += aboutbox_test.h
