include(../../variables.pri)

TEMPLATE = lib
LANGUAGE = C
TARGET   = btrack

#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

#CONFIG      += plugin

INCLUDEPATH += ../samplerate/include
DEPENDPATH += ../samplerate

LIBS        += ../samplerate/libsamplerate.so

DEFINES += USE_FFTW

DISTFILES += \
    LICENSE.txt \
    README.md

HEADERS += \
    src/BTrack.h \
    src/CircularBuffer.h \
    src/OnsetDetectionFunction.h

SOURCES += \
    src/BTrack.cpp \
    src/OnsetDetectionFunction.cpp

