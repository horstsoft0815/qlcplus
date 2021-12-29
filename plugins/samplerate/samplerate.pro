include(../../variables.pri)

TEMPLATE = lib
LANGUAGE = C
TARGET   = samplerate

#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

#CONFIG      += plugin

target.path = $$INSTALLROOT/$$PLUGINDIR
INSTALLS   += target


DISTFILES += \
    AUTHORS \
    COPYING \
    README.md

HEADERS += \
    config.h \
    include/samplerate.h \
    src/common.h \
    src/fastest_coeffs.h \
    src/high_qual_coeffs.h \
    src/mid_qual_coeffs.h

SOURCES += \
    src/samplerate.c \
    src/src_linear.c \
    src/src_sinc.c \
    src/src_zoh.c
