include(../../variables.pri)

TEMPLATE = lib
LANGUAGE = C++
TARGET   = btrack_vamp

#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

#CONFIG      += plugin

DEFINES += USE_FFTW

DISTFILES += \
    INSTALL.md \
    vamp-plugin.list \
    vamp-plugin.map

HEADERS += \
    BTrackVamp.h

SOURCES += \
    BTrackVamp.cpp \
    plugins.cpp



