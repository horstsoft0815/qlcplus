include(../../variables.pri)

TEMPLATE = lib
LANGUAGE = C++
TARGET   = analyzer
QT += core gui
QT += qml
QT += xml
QT += widgets
QT += svg

CONFIG += c++17
CONFIG += plugin

#QTPLUGIN =

INCLUDEPATH += ../../plugins/qm-dsp
DEPENDPATH += ../../plugins/qm-dsp

LIBS        += ../../plugins/qm-dsp/libqm-dsp.so

HEADERS += \
    analyzer/constants.h \
    analyzer/plugins/analyzerplugin.h \
    analyzer/plugins/analyzerqueenmarybeats.h \
    analyzer/plugins/buffering_utils.h \
    audio/frame.h \
    audio/signalinfo.h \
    audio/types.h \
    engine/engine.h \
    proto/beats.pb.h \
    track/beats.h \
    track/beatutils.h \
    track/bpm.h \
    util/assert.h \
    util/fpclassify.h \
    util/macros.h \
    util/math.h \
    util/memory.h \
    util/optional.h \
    util/platform.h \
    util/sample.h \
    util/sample_autogen.h \
    util/samplebuffer.h \
    util/types.h

SOURCES += \
    analyzer/plugins/analyzerqueenmarybeats.cpp \
    analyzer/plugins/buffering_utils.cpp \
    audio/frame.cpp \
    audio/signalinfo.cpp \
    audio/types.cpp \
    proto/beats.pb.cc \
    track/beats.cpp \
    track/beatutils.cpp \
    track/bpm.cpp \
    util/fpclassify.cpp \
    util/sample.cpp \
    util/samplebuffer.cpp


