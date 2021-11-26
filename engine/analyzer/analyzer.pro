include(../../variables.pri)

TEMPLATE = lib
LANGUAGE = C++
TARGET   = analyzer

QT += core gui
QT += qml
QT += xml
QT += widgets
QT += svg


CONFIG += plugin

#QTPLUGIN =

INCLUDEPATH += ../../plugins/qm-dsp
DEPENDPATH += ../../plugins/qm-dsp

LIBS        += -L../../plugins/qm-dsp -lqm-dsp.so


