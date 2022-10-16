include(../../../variables.pri)

TEMPLATE = lib
LANGUAGE = C++
TARGET   = qlcplusaudio
CONFIG  += staticlib

QT      += core
greaterThan(QT_MAJOR_VERSION, 4) {
  QT += multimedia
  macx:QT_CONFIG -= no-pkg-config
  win32:QT += widgets
}

CONFIG += link_pkgconfig

INCLUDEPATH += ../../src ../../../plugins/interfaces

# Beat tracking
INCLUDEPATH += ../plugins/kissfft/
INCLUDEPATH += ../plugins/samplerate/include

DEPENDPATH  += ../plugins/kissfft
DEPENDPATH  += ../plugins/samplerate
QMAKE_LIBDIR+= ../plugins/kissfft
QMAKE_LIBDIR+= ../plugins/samplerate

LIBS        +=  -lkissfft -lsamplerate
DEFINES += USE_KISS_FFT

HEADERS += audio.h \
           BTrack.h \
           CircularBuffer.h \
           OnsetDetectionFunction.h \
           audiodecoder.h \
           audiorenderer.h \
           audioparameters.h \
           audiocapture.h \
           audioplugincache.h

lessThan(QT_MAJOR_VERSION, 5) {
  unix:!macx:HEADERS += audiorenderer_alsa.h audiocapture_alsa.h
  win32:HEADERS += audiorenderer_waveout.h audiocapture_wavein.h
}
else {
  HEADERS += audiorenderer_qt.h audiocapture_qt.h
}

SOURCES += audio.cpp \
           BTrack.cpp \
           CircularBuffer.cpp \
           OnsetDetectionFunction.cpp \
           audiodecoder.cpp \
           audiorenderer.cpp \
           audioparameters.cpp \
           audiocapture.cpp \
           audioplugincache.cpp

lessThan(QT_MAJOR_VERSION, 5) {
  unix:!macx:SOURCES += audiorenderer_alsa.cpp audiocapture_alsa.cpp
  win32:SOURCES += audiorenderer_waveout.cpp audiocapture_wavein.cpp

  macx {
    system(pkg-config --exists portaudio-2.0) {
      DEFINES += HAS_PORTAUDIO
      PKGCONFIG += portaudio-2.0
      HEADERS += audiorenderer_portaudio.h audiocapture_portaudio.h
      SOURCES += audiorenderer_portaudio.cpp audiocapture_portaudio.cpp
    }
  }
}
else {
  SOURCES += audiorenderer_qt.cpp audiocapture_qt.cpp
}

!android:!ios {
  system(pkg-config --exists fftw3) {
    DEFINES += HAS_FFTW3
    PKGCONFIG += fftw3
    macx:LIBS += -lfftw3
  }

  unix:!macx:LIBS += -lasound
}
