include(../../../../variables.pri)

TEMPLATE = lib
LANGUAGE = C++
TARGET   = kissfft

INCLUDEPATH += ../../src
CONFIG  += staticlib

CONFIG += link_pkgconfig
#PKGCONFIG   += kissfft

#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

#CONFIG      += plugin

target.path = $$INSTALLROOT/$$PLUGINDIR
INSTALLS   += target

DISTFILES += \
    COPYING \
    LICENSES/BSD-3-Clause \
    LICENSES/Unlicense \
    README.md \
    README.simd \
    TIPS

SOURCES += \
    kfc.c \
    kiss_fft.c \
    kiss_fftnd.c \
    kiss_fftndr.c \
    kiss_fftr.c

HEADERS += \
    _kiss_fft_guts.h \
    kfc.h \
    kiss_fft.h \
    kiss_fft_log.h \
    kiss_fftnd.h \
    kiss_fftndr.h \
    kiss_fftr.h \
    kissfft.hh \
    kissfft_i32.hh
