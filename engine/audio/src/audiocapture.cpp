/*
  Q Light Controller Plus
  audiocapture.cpp

  Copyright (c) Massimo Callegari
  based on libbeat code by Maximilian Güntner

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

#include <QSettings>
#include <QDebug>
#include <qmath.h>
#include <iostream>

#include "audiocapture.h"
#include "BTrack.h"

#ifdef HAS_FFTW3
#include "fftw3.h"
#endif

#define USE_HANNING
#define CLEAR_FFT_NOISE
#define M_2PI       6.28318530718           /* 2*pi */

#include "kiss_fft.h"

#ifdef _0
enum OnsetDetectionFunctionTypeTest
{
    eEnergyEnvelope,
    eEnergyDifference,
    eSpectralDifference,
    eSpectralDifferenceHWR,
    ePhaseDeviation,
    eComplexSpectralDifference,
    eComplexSpectralDifferenceHWR,
    eHighFrequencyContent,
    eHighFrequencySpectralDifference,
    eHighFrequencySpectralDifferenceHWR
};

enum WindowTypeTest
{
    eRectangularWindow,
    eHanningWindow,
    eHammingWindow,
    eBlackmanWindow,
    eTukeyWindow
};

class OnsetDetectionFunctionTest
{
public:
    OnsetDetectionFunctionTest (int hopSize_p, int frameSize_p);

    OnsetDetectionFunctionTest (int hopSize_p, int frameSize_p, int onsetDetectionFunctionType_p, int windowType_p);

    /** Destructor */
    ~OnsetDetectionFunctionTest();

    /** Initialisation function for only updating hop size and frame size (and not window type
     * or onset detection function type
     * @param hopSize_p the hop size in audio samples
     * @param frameSize_p the frame size in audio samples
     */
    void initialise (int hopSize_p, int frameSize_p);

    /** Initialisation Function
     * @param hopSize_p the hop size in audio samples
     * @param frameSize_p the frame size in audio samples
     * @param onsetDetectionFunctionType_p the type of onset detection function to use - (see OnsetDetectionFunctionType)
     * @param windowType_p the type of window to use (see WindowType)
     */
    void initialise (int hopSize_p, int frameSize_p, int onsetDetectionFunctionType_p, int windowType_p);
    void initialiseFFT();
    void freeFFT();

    /** Calculate a Rectangular window */
    void calculateRectangularWindow();

    /** Calculate a Hanning window */
    void calculateHanningWindow();

    /** Calculate a Hamming window */
    void calclulateHammingWindow();

    /** Calculate a Blackman window */
    void calculateBlackmanWindow();

    /** Calculate a Tukey window */
    void calculateTukeyWindow();

    double m_Pi;							/**< pi, the constant */

    int m_FrameSize;						/**< audio framesize */
    int m_HopSize;						/**< audio hopsize */
    int m_OnsetDetectionFunctionType;		/**< type of detection function */
    int m_WindowType;                     /**< type of window used in calculations */



    kiss_fft_cfg m_Cfg;                   /**< Kiss FFT configuration */
    kiss_fft_cpx* m_FFTIn;                /**< FFT input samples, in complex form */
    kiss_fft_cpx* m_FFTOut;               /**< FFT output samples, in complex form */
    std::vector<std::vector<double> > m_ComplexOut;


    //=======================================================================
    bool m_Initialised = false;					/**< flag indicating whether buffers and FFT plans are initialised */

    std::vector<double> m_Frame;          /**< audio frame */
    std::vector<double> m_Window;         /**< window */

    double m_PrevEnergySum;				/**< to hold the previous energy sum value */

    std::vector<double> m_MagSpec;        /**< magnitude spectrum */
    std::vector<double> m_PrevMagSpec;    /**< previous magnitude spectrum */

    std::vector<double> m_Phase;          /**< FFT phase values */
    std::vector<double> m_PrevPhase;      /**< previous phase values */
    std::vector<double> m_PrevPhase2;     /**< second order previous phase values */

};

OnsetDetectionFunctionTest::OnsetDetectionFunctionTest (
        int hopSize_p,
        int frameSize_p)
 :  m_OnsetDetectionFunctionType (eComplexSpectralDifferenceHWR),
   m_WindowType (eHanningWindow)
{
    // indicate that we have not initialised yet
    m_Initialised = false;

    // set pi
    m_Pi = 3.14159265358979;

    // initialise with arguments to constructor
    initialise (hopSize_p, frameSize_p, eComplexSpectralDifferenceHWR, eHanningWindow);
}

//=======================================================================
OnsetDetectionFunctionTest::OnsetDetectionFunctionTest(
        int hopSize_p,
        int frameSize_p,
        int onsetDetectionFunctionType_p,
        int windowType_p)
 :  m_OnsetDetectionFunctionType (eComplexSpectralDifferenceHWR), m_WindowType (eHanningWindow)
{
    // indicate that we have not initialised yet
    m_Initialised = false;

    // set pi
    m_Pi = 3.14159265358979;

    // initialise with arguments to constructor
    initialise (hopSize_p, frameSize_p, onsetDetectionFunctionType_p, windowType_p);
}


//=======================================================================
OnsetDetectionFunctionTest::~OnsetDetectionFunctionTest()
{
    if (m_Initialised)
    {
        //freeFFT();
    }
}

//=======================================================================
void OnsetDetectionFunctionTest::initialise (int hopSize_p, int frameSize_p)
{
    // use the already initialised onset detection function and window type and
    // pass the new frame and hop size to the main initialisation function
    initialise (hopSize_p, frameSize_p, m_OnsetDetectionFunctionType, m_WindowType);
}

//=======================================================================
void OnsetDetectionFunctionTest::initialise(
        int hopSize_p,
        int frameSize_p,
        int onsetDetectionFunctionType_p,
        int windowType_p)
{
    m_HopSize = hopSize_p; // set hopsize
    m_FrameSize = frameSize_p; // set framesize

    m_OnsetDetectionFunctionType = onsetDetectionFunctionType_p; // set detection function type
    m_WindowType = windowType_p; // set window type

    // initialise buffers
    m_Frame.resize (m_FrameSize);
    m_Window.resize (m_FrameSize);
    m_MagSpec.resize (m_FrameSize);
    m_PrevMagSpec.resize (m_FrameSize);
    m_Phase.resize (m_FrameSize);
    m_PrevPhase.resize (m_FrameSize);
    m_PrevPhase2.resize (m_FrameSize);


    // set the window to the specified type
    switch (m_WindowType)
    {
        case eRectangularWindow:
            calculateRectangularWindow();		// Rectangular window
            break;
        case eHanningWindow:
            calculateHanningWindow();			// Hanning Window
            break;
        case eHammingWindow:
            calclulateHammingWindow();			// Hamming Window
            break;
        case eBlackmanWindow:
            calculateBlackmanWindow();			// Blackman Window
            break;
        case eTukeyWindow:
            calculateTukeyWindow();             // Tukey Window
            break;
        default:
            calculateHanningWindow();			// DEFAULT: Hanning Window
    }

    // initialise previous magnitude spectrum to zero
    for (int i = 0; i < m_FrameSize; i++)
    {
        m_PrevMagSpec[i] = 0.0;
        m_PrevPhase[i] = 0.0;
        m_PrevPhase2[i] = 0.0;
        m_Frame[i] = 0.0;
    }

    m_PrevEnergySum = 0.0;	// initialise previous energy sum value to zero

    initialiseFFT();
}

//=======================================================================
void OnsetDetectionFunctionTest::initialiseFFT()
{
    if (m_Initialised) // if we have already initialised FFT plan
    {
        freeFFT();
    }

#ifdef USE_FFTW
    m_ComplexIn = (fftw_complex*) fftw_malloc (sizeof(fftw_complex) * m_FrameSize);		// complex array to hold fft data
    m_ComplexOut = (fftw_complex*) fftw_malloc (sizeof(fftw_complex) * m_FrameSize);	// complex array to hold fft data
    //m_Plan = fftw_plan_dft_1d (m_FrameSize, m_ComplexIn, m_ComplexOut, FFTW_FORWARD, FFTW_ESTIMATE);	// FFT plan initialisation
#endif

#ifdef USE_KISS_FFT
    m_ComplexOut.resize (m_FrameSize);

    for (int i = 0; i < m_FrameSize;i++)
    {
        m_ComplexOut[i].resize(2);
    }

    m_FFTIn = new kiss_fft_cpx[m_FrameSize];
    m_FFTOut = new kiss_fft_cpx[m_FrameSize];
    m_Cfg = kiss_fft_alloc (m_FrameSize, 0, 0, 0);
#endif

    m_Initialised = true;
}

//=======================================================================
void OnsetDetectionFunctionTest::freeFFT()
{
#ifdef USE_FFTW
    //fftw_destroy_plan (m_Plan);
    fftw_free (m_ComplexIn);
    fftw_free (m_ComplexOut);
#endif

#ifdef USE_KISS_FFT
    free (m_Cfg);
    delete [] m_FFTIn;
    delete [] m_FFTOut;
#endif
}


void OnsetDetectionFunctionTest::calculateHanningWindow()
{
    double N;		// variable to store framesize minus 1

    N = (double) (m_FrameSize-1);	// framesize minus 1

    // Hanning window calculation
    for (int n = 0; n < m_FrameSize; n++)
    {
        m_Window[n] = 0.5 * (1 - cos (2 * m_Pi * (n / N)));
    }
}

//=======================================================================
void OnsetDetectionFunctionTest::calclulateHammingWindow()
{
    double N;		// variable to store framesize minus 1
    double n_val;	// double version of index 'n'

    N = (double) (m_FrameSize-1);	// framesize minus 1
    n_val = 0;

    // Hamming window calculation
    for (int n = 0;n < m_FrameSize;n++)
    {
        m_Window[n] = 0.54 - (0.46 * cos (2 * m_Pi * (n_val/N)));
        n_val = n_val+1;
    }
}

//=======================================================================
void OnsetDetectionFunctionTest::calculateBlackmanWindow()
{
    double N;		// variable to store framesize minus 1
    double n_val;	// double version of index 'n'

    N = (double) (m_FrameSize-1);	// framesize minus 1
    n_val = 0;

    // Blackman window calculation
    for (int n = 0;n < m_FrameSize;n++)
    {
        m_Window[n] = 0.42 - (0.5*cos(2*m_Pi*(n_val/N))) + (0.08*cos(4*m_Pi*(n_val/N)));
        n_val = n_val+1;
    }
}

//=======================================================================
void OnsetDetectionFunctionTest::calculateTukeyWindow()
{
    double N;		// variable to store framesize minus 1
    double n_val;	// double version of index 'n'
    double alpha;	// alpha [default value = 0.5];

    alpha = 0.5;

    N = (double) (m_FrameSize-1);	// framesize minus 1

    // Tukey window calculation

    n_val = (double) (-1*((m_FrameSize/2)))+1;

    for (int n = 0;n < m_FrameSize;n++)	// left taper
    {
        if ((n_val >= 0) && (n_val <= (alpha*(N/2))))
        {
            m_Window[n] = 1.0;
        }
        else if ((n_val <= 0) && (n_val >= (-1*alpha*(N/2))))
        {
            m_Window[n] = 1.0;
        }
        else
        {
            m_Window[n] = 0.5*(1+cos(m_Pi*(((2*n_val)/(alpha*N))-1)));
        }

        n_val = n_val+1;
    }

}

//=======================================================================
void OnsetDetectionFunctionTest::calculateRectangularWindow()
{
    // Rectangular window calculation
    for (int n = 0;n < m_FrameSize;n++)
    {
        m_Window[n] = 1.0;
    }
}
#endif



AudioCapture::AudioCapture (QObject* parent)
    : QThread (parent)
    , m_IsFrequencyAnalysisActive(true)
    , m_IsBeatAnalysisActive(false)
    , m_userStop(true)
    , m_pause(false)
    , m_bufferSize(AUDIO_DEFAULT_BUFFER_SIZE)
    , m_captureSize(0)
    , m_sampleRate(AUDIO_DEFAULT_SAMPLE_RATE)
    , m_channels(AUDIO_DEFAULT_CHANNELS)
    , m_audioBuffer(NULL)
    , m_audioMixdown(NULL)
    , m_fftInputBuffer(m_bufferSize)
    , m_fftOutputBuffer(NULL)
    , m_AudioBufferF64(m_bufferSize)
    , m_PrevAudioBufferF64()
    , m_FrameSize(m_bufferSize)
    , m_HopSize(m_FrameSize/2)
    , m_BTrack(QSharedPointer<BTrack>(new BTrack(2*m_HopSize, 2*m_FrameSize)))
    //, m_BTrack(m_HopSize, m_FrameSize)
{
    QSettings settings;
    QVariant var = settings.value(SETTINGS_AUDIO_INPUT_SRATE);

    if (var.isValid() == true)
        m_sampleRate = var.toInt();

    var = settings.value(SETTINGS_AUDIO_INPUT_CHANNELS);

    if (var.isValid() == true)
        m_channels = var.toInt();

    qDebug() << "[AudioCapture] initialize" << m_sampleRate << m_channels;

    m_captureSize = m_bufferSize * m_channels;

    m_audioBuffer = new int16_t[m_captureSize];
    m_audioMixdown = new int16_t[m_bufferSize];
#ifdef HAS_FFTW3
    m_fftOutputBuffer = fftw_malloc(sizeof(fftw_complex) * m_bufferSize);
#endif
}

AudioCapture::~AudioCapture()
{
    // stop() has to be called from the implementation class
    Q_ASSERT(!this->isRunning());

    delete[] m_audioBuffer;
    delete[] m_audioMixdown;
#ifdef HAS_FFTW3
    if (m_fftOutputBuffer)
        fftw_free(m_fftOutputBuffer);
#endif
}

int AudioCapture::defaultBarsNumber()
{
    return FREQ_SUBBANDS_DEFAULT_NUMBER;
}

void AudioCapture::registerBandsNumber(int number)
{
    qDebug() << "[AudioCapture] registering" << number << "bands";

    QMutexLocker locker(&m_mutex);

    bool firstBand = m_fftMagnitudeMap.isEmpty();
    if (number > 0 && number <= FREQ_SUBBANDS_MAX_NUMBER)
    {
        if (m_fftMagnitudeMap.contains(number) == false)
        {
            BandsData newBands;
            newBands.m_registerCounter = 1;
            newBands.m_fftMagnitudeBuffer = QVector<double>(number);
            m_fftMagnitudeMap[number] = newBands;
        }
        else
            m_fftMagnitudeMap[number].m_registerCounter++;

        if (firstBand)
        {
            locker.unlock();
            start();
        }
    }
}

void AudioCapture::unregisterBandsNumber(int number)
{
    qDebug() << "[AudioCapture] unregistering" << number << "bands";

    QMutexLocker locker(&m_mutex);

    if (m_fftMagnitudeMap.contains(number))
    {
        m_fftMagnitudeMap[number].m_registerCounter--;
        if (m_fftMagnitudeMap[number].m_registerCounter == 0)
            m_fftMagnitudeMap.remove(number);

        if (m_fftMagnitudeMap.isEmpty())
        {
            locker.unlock();
            stop();
        }
    }
}

void AudioCapture::stop()
{
    qDebug() << "[AudioCapture] stop capture";
    while (this->isRunning())
    {
        m_userStop = true;
        usleep(10000);
    }
}

double AudioCapture::fillBandsData(int number)
{
    // m_fftOutputBuffer contains the real and imaginary data of a spectrum
    // representing all the frequencies from 0 to m_sampleRate Hz.
    // I will just consider 0 to 5000Hz and will calculate average magnitude
    // for the number of desired bands.
    double maxMagnitude = 0.;
#ifdef HAS_FFTW3
    unsigned int i = 1; // skip DC bin
    int subBandWidth = ((m_bufferSize * SPECTRUM_MAX_FREQUENCY) / m_sampleRate) / number;

    for (int b = 0; b < number; b++)
    {
        double magnitudeSum = 0.;
        for (int s = 0; s < subBandWidth; s++, i++)
        {
            if (i == m_bufferSize)
                break;
            magnitudeSum += qSqrt((((fftw_complex*)m_fftOutputBuffer)[i][0] * ((fftw_complex*)m_fftOutputBuffer)[i][0]) +
                                  (((fftw_complex*)m_fftOutputBuffer)[i][1] * ((fftw_complex*)m_fftOutputBuffer)[i][1]));
        }
        double bandMagnitude = (magnitudeSum / (subBandWidth * M_2PI));
        m_fftMagnitudeMap[number].m_fftMagnitudeBuffer[b] = bandMagnitude;
        if (maxMagnitude < bandMagnitude)
            maxMagnitude = bandMagnitude;
    }
#else
    Q_UNUSED(number)
#endif
    return maxMagnitude;
}

void AudioCapture::processData()
{
#ifdef HAS_FFTW3
    unsigned int i, j;
    double pwrSum = 0.;
    double maxMagnitude = 0.;

    // 1 ********* Initialize FFTW
    fftw_plan plan_forward;
    plan_forward = fftw_plan_dft_r2c_1d(m_bufferSize, m_fftInputBuffer.data(), (fftw_complex*)m_fftOutputBuffer , 0);

    // 2 ********* Apply a window to audio data
    // *********** and convert it to doubles

    // Mix down the channels to mono
    for (i = 0; i < m_bufferSize; i++)
    {
        m_audioMixdown[i] = 0;
        for (j = 0; j < m_channels; j++)
        {
            m_audioMixdown[i] += m_audioBuffer[i*m_channels + j] / m_channels;
        }
    }

    for (i = 0; i < m_bufferSize; i++)
    {
#ifdef USE_BLACKMAN
        double a0 = (1-0.16)/2;
        double a1 = 0.5;
        double a2 = 0.16/2;
        m_fftInputBuffer[i] = m_audioMixdown[i] * (a0 - a1 * qCos((M_2PI * i) / (m_bufferSize - 1)) +
                              a2 * qCos((2 * M_2PI * i) / (m_bufferSize - 1))) / 32768.;
#endif
#ifdef USE_HANNING
        m_fftInputBuffer[i] = m_audioMixdown[i] * (0.5 * (1.00 - qCos((M_2PI * i) / (m_bufferSize - 1)))) / 32768.;
#endif
#ifdef USE_NO_WINDOW
        m_fftInputBuffer[i] = (double)m_audioMixdown[i] / 32768.;
#endif
    }

    // 3 ********* Perform FFT
    fftw_execute(plan_forward);
    fftw_destroy_plan(plan_forward);

    // 4 ********* Clear FFT noise
#ifdef CLEAR_FFT_NOISE
    //We delete some values since these will ruin our output
    for (unsigned int n = 0; n < 5 && n < m_bufferSize; n++)
    {
        ((fftw_complex*)m_fftOutputBuffer)[n][0] = 0;
        ((fftw_complex*)m_fftOutputBuffer)[n][1] = 0;
    }
#endif

    // 5 ********* Calculate the average signal power
    foreach(int barsNumber, m_fftMagnitudeMap.keys())
    {
        maxMagnitude = fillBandsData(barsNumber);
        pwrSum = 0.;
        for (int n = 0; n < barsNumber; n++)
        {
            pwrSum += m_fftMagnitudeMap[barsNumber].m_fftMagnitudeBuffer[n];
        }
        m_signalPower = 32768 * pwrSum * qSqrt(M_2PI) / (double)barsNumber;
        emit dataProcessed(m_fftMagnitudeMap[barsNumber].m_fftMagnitudeBuffer.data(),
                           m_fftMagnitudeMap[barsNumber].m_fftMagnitudeBuffer.size(),
                           maxMagnitude, m_signalPower);
    }
#endif
}

#if 0
void AudioCapture::detectBeat(QSharedPointer<BTrack> bTrack_p)
{
    if(bTrack_p)
    {
        assert(m_bufferSize == m_AudioBufferF64.size());

        // Mix down the channels to mono
        for (unsigned int i = 0; i < m_bufferSize; i++)
        {
            m_AudioBufferF64[i] = 0;
            for (unsigned int j = 0; j < m_channels; j++)
            {
                m_AudioBufferF64[i] += static_cast<double>(m_audioBuffer[i*m_channels + j]) / m_channels;
            }
        }

        if(m_PrevAudioBufferF64.empty())
        {
            m_PrevAudioBufferF64 = m_AudioBufferF64;
            bTrack_p->processAudioFrame(m_AudioBufferF64.data());

            emit detectedBeat(m_BTrack->beatDueInCurrentFrame());
        }
        else
        {
            bool isBeat{false};

            std::vector<double> combinedFrame(m_PrevAudioBufferF64);
            std::copy(m_AudioBufferF64.begin(), m_AudioBufferF64.end(), std::back_inserter(combinedFrame));

            std::vector<double> frameWithOverlap(m_FrameSize);
            std::size_t startPosition = m_HopSize;
            double* pStart = combinedFrame.data() + startPosition;
            double* pEnd = pStart + m_FrameSize;
            while(pEnd < (combinedFrame.data() + combinedFrame.size()) && m_HopSize>0)
            {
                std::copy(pStart, pEnd, frameWithOverlap.begin());

                bTrack_p->processAudioFrame(frameWithOverlap.data());

                if(bTrack_p->beatDueInCurrentFrame())
                {
                    isBeat = true;
                    break;
                }

                startPosition += m_HopSize;
                pStart = combinedFrame.data() + startPosition;
                pEnd = pStart + m_FrameSize;
            }

            if(isBeat)
            {
                std::cout << "beat\n";
            }

            emit detectedBeat(isBeat);
        }
    } // if(bTrack_p)
}
#endif

void AudioCapture::detectBeat(QSharedPointer<BTrack> bTrack_p)
{
    if(bTrack_p)
    {
        assert(m_bufferSize == m_AudioBufferF64.size());

        // Mix down the channels to mono
        for (unsigned int i = 0; i < m_bufferSize; i++)
        {
            m_AudioBufferF64[i] = 0;
            for (unsigned int j = 0; j < m_channels; j++)
            {
                m_AudioBufferF64[i] += static_cast<double>(m_audioBuffer[i*m_channels + j]) / m_channels;
            }
        }

        bool isBeat{false};
        bTrack_p->processAudioFrame(m_AudioBufferF64.data());

        if(bTrack_p->beatDueInCurrentFrame())
        {
            isBeat = true;
        }

        const quint32 tempo = static_cast<quint32>(bTrack_p->getTempo()); // limited < 160, hence no overflow problems

        emit detectedBeat(isBeat);
        emit detectedBPM(tempo);
    } // if(bTrack_p)
}



void AudioCapture::run()
{
    qDebug() << "[AudioCapture] start capture";

    m_userStop = false;

    if (!initialize())
    {
        qWarning() << "[AudioCapture] Could not initialize audio capture, abandon";
        return;
    }

    //BTrack m_BTrack(m_HopSize, m_FrameSize);
    //OnsetDetectionFunction m_Odf(m_HopSize, m_FrameSize);
    //kiss_fft_cfg m_Cfg = kiss_fft_alloc (m_FrameSize, 0, 0, 0);
    //kiss_fft_cpx* m_FFTIn = new kiss_fft_cpx[m_FrameSize];

    while (!m_userStop)
    {
        if (m_pause == false && m_captureSize != 0)
        {
            if (readAudio(m_captureSize) == true)
            {
                QMutexLocker locker(&m_mutex);

                if(m_IsFrequencyAnalysisActive)
                {
                    processData();
                }

                if(m_IsBeatAnalysisActive)
                {

                    detectBeat(m_BTrack);
                }
            }
            else
            {
                //qDebug() << "Error reading data from audio source";
                QThread::msleep(5);
            }

        }
        else
        {
            QThread::msleep(15);
        }

        QThread::yieldCurrentThread();
    }

    //std::cout << m_Cfg;
    //std::cout << m_FFTIn;

    uninitialize();
}
