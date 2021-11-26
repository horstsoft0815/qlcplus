include(../../variables.pri)

TEMPLATE = lib
LANGUAGE = C++
TARGET   = qm-dsp

#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG      += plugin

DEFINES += "kiss_fft_scalar=double"

# This must be after "TARGET = " and before target installation so that
# install_name_tool can be run before target installation
#macx:include(../../platforms/macos/nametool.pri)

target.path = $$INSTALLROOT/$$PLUGINDIR
INSTALLS   += target

INCLUDEPATH += base
INCLUDEPATH += dsp
INCLUDEPATH += ext
INCLUDEPATH += hmm
INCLUDEPATH += include
INCLUDEPATH += maths
               
INCLUDEPATH += dsp\chromagram
INCLUDEPATH += dsp\keydetection
INCLUDEPATH += dsp\mfcc
INCLUDEPATH += dsp\onsets
INCLUDEPATH += dsp\phasevocoder
INCLUDEPATH += dsp\rateconversion
INCLUDEPATH += dsp\rhythm
INCLUDEPATH += dsp\segmentation
INCLUDEPATH += dsp\signalconditioning
INCLUDEPATH += dsp\tempotracking
INCLUDEPATH += dsp\tonal
INCLUDEPATH += dsp\transforms
INCLUDEPATH += dps\wavelet
               
INCLUDEPATH += ext\kissfft
INCLUDEPATH += ext\kissfft\tools
               
INCLUDEPATH += maths\pca

#win32:LIBS += -lsetupapi -lwinmm -lhid
#macx:LIBS += -framework IOKit -framework CoreFoundation

#HEADERS += ../interfaces/qlcioplugin.h

DISTFILES += \
    CONTRIBUTING.md \
    COPYING \
    README.md \
    ext/kissfft/CHANGELOG \
    ext/kissfft/COPYING \
    ext/kissfft/README \
    ext/kissfft/README.simd \
    ext/kissfft/TIPS \
    maths/pca/data.txt \
    mixxx-changes.patch

HEADERS += \
    base/KaiserWindow.h \
    base/Pitch.h \
    base/Restrict.h \
    base/SincWindow.h \
    base/Window.h \
    dsp/chromagram/Chromagram.h \
    dsp/chromagram/ConstantQ.h \
    dsp/keydetection/GetKeyMode.h \
    dsp/mfcc/MFCC.h \
    dsp/onsets/DetectionFunction.h \
    dsp/onsets/PeakPicking.h \
    dsp/phasevocoder/PhaseVocoder.h \
    dsp/rateconversion/Decimator.h \
    dsp/rateconversion/DecimatorB.h \
    dsp/rateconversion/Resampler.h \
    dsp/rhythm/BeatSpectrum.h \
    dsp/segmentation/ClusterMeltSegmenter.h \
    dsp/segmentation/Segmenter.h \
    dsp/segmentation/cluster_melt.h \
    dsp/segmentation/cluster_segmenter.h \
    dsp/segmentation/segment.h \
    dsp/signalconditioning/DFProcess.h \
    dsp/signalconditioning/FiltFilt.h \
    dsp/signalconditioning/Filter.h \
    dsp/signalconditioning/Framer.h \
    dsp/tempotracking/DownBeat.h \
    dsp/tempotracking/TempoTrack.h \
    dsp/tempotracking/TempoTrackV2.h \
    dsp/tonal/ChangeDetectionFunction.h \
    dsp/tonal/TCSgram.h \
    dsp/tonal/TonalEstimator.h \
    dsp/transforms/DCT.h \
    dsp/transforms/FFT.h \
    dsp/wavelet/Wavelet.h \
    ext/kissfft/_kiss_fft_guts.h \
    ext/kissfft/kiss_fft.h \
    ext/kissfft/kissfft.hh \
    ext/kissfft/tools/kiss_fftr.h \
    hmm/hmm.h \
    include/cblas.h \
    include/clapack.h \
    maths/Correlation.h \
    maths/CosineDistance.h \
    maths/KLDivergence.h \
    maths/MathAliases.h \
    maths/MathUtilities.h \
    maths/MedianFilter.h \
    maths/Polyfit.h \
    maths/nan-inf.h \
    maths/pca/pca.h

SOURCES += \
    base/KaiserWindow.cpp \
    base/Pitch.cpp \
    base/SincWindow.cpp \
    dsp/chromagram/Chromagram.cpp \
    dsp/chromagram/ConstantQ.cpp \
    dsp/keydetection/GetKeyMode.cpp \
    dsp/mfcc/MFCC.cpp \
    dsp/onsets/DetectionFunction.cpp \
    dsp/onsets/PeakPicking.cpp \
    dsp/phasevocoder/PhaseVocoder.cpp \
    dsp/rateconversion/Decimator.cpp \
    dsp/rateconversion/DecimatorB.cpp \
    dsp/rateconversion/Resampler.cpp \
    dsp/rhythm/BeatSpectrum.cpp \
    dsp/segmentation/ClusterMeltSegmenter.cpp \
    dsp/segmentation/Segmenter.cpp \
    dsp/segmentation/cluster_melt.c \
    dsp/segmentation/cluster_segmenter.c \
    dsp/signalconditioning/DFProcess.cpp \
    dsp/signalconditioning/FiltFilt.cpp \
    dsp/signalconditioning/Filter.cpp \
    dsp/signalconditioning/Framer.cpp \
    dsp/tempotracking/DownBeat.cpp \
    dsp/tempotracking/TempoTrack.cpp \
    dsp/tempotracking/TempoTrackV2.cpp \
    dsp/tonal/ChangeDetectionFunction.cpp \
    dsp/tonal/TCSgram.cpp \
    dsp/tonal/TonalEstimator.cpp \
    dsp/transforms/DCT.cpp \
    dsp/transforms/FFT.cpp \
    dsp/wavelet/Wavelet.cpp \
    ext/kissfft/kiss_fft.c \
    ext/kissfft/tools/kiss_fftr.c \
    hmm/hmm.c \
    maths/Correlation.cpp \
    maths/CosineDistance.cpp \
    maths/KLDivergence.cpp \
    maths/MathUtilities.cpp \
    maths/pca/pca.c

