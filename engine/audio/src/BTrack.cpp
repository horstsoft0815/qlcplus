//=======================================================================
/** @file BTrack.cpp
 *  @brief BTrack - a real-time beat tracker
 *  @author Adam Stark
 *  @copyright Copyright (C) 2008-2014  Queen Mary University of London
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//=======================================================================

#include "BTrack.h"
#include "samplerate.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>


#define CIRCULARBUFFERFAC 1
#define MINBPM 80
#define MAXBPM 180

// firstly make sure tempo is between min and max bpm..
static void ReduceTempo(double& tempo_p)
{
    while (tempo_p > MAXBPM)
    {
        tempo_p = tempo_p/2;
    }
#if 0
    while (tempo_p < MINBPM)
    {
        tempo_p = tempo_p * 2;
    }
#endif
}

//=======================================================================
BTrack::BTrack()
 :  m_Odf (512, 1024, eComplexSpectralDifferenceHWR, eHanningWindow)
{
    initialise (512, 1024);
}

//=======================================================================
BTrack::BTrack (int hopSize_p)
 :  m_Odf(hopSize_p, 2*hopSize_p, eComplexSpectralDifferenceHWR, eHanningWindow)
{	
    initialise (hopSize_p, 2*hopSize_p);
}

//=======================================================================
BTrack::BTrack (int hopSize_p, int frameSize_p)
 : m_Odf (hopSize_p, frameSize_p, eComplexSpectralDifferenceHWR, eHanningWindow)
{
    initialise (hopSize_p, frameSize_p);
}

//=======================================================================
BTrack::~BTrack()
{
#ifdef USE_FFTW
    // destroy fft plan
    //fftw_destroy_plan (acfForwardFFT);
    //fftw_destroy_plan (acfBackwardFFT);
    fftw_free (m_ComplexIn);
    fftw_free (m_ComplexOut);
#endif
    
#ifdef USE_KISS_FFT
    /*free (m_CfgForwards);
    free (m_CfgBackwards);
    delete [] m_FFTIn;
    delete [] m_FFTOut;*/
#endif
}

//=======================================================================
double BTrack::getBeatTimeInSeconds (long frameNumber_p, int hopSize_p, int fs_p)
{
    double hop = (double) hopSize_p;
    double samplingFrequency = (double) fs_p;
    double frameNum = (double) frameNumber_p;
    
    return ((hop / samplingFrequency) * frameNum);
}

//=======================================================================
double BTrack::getBeatTimeInSeconds (int frameNumber_p, int hopSize_p, int fs_p)
{
    long frameNum = (long) frameNumber_p;
    
    return getBeatTimeInSeconds (frameNum, hopSize_p, fs_p);
}



//=======================================================================
void BTrack::initialise (int hopSize_p, int frameSize_p)
{
    double rayparam = 43;
	double pi = 3.14159265;
	
	
	// initialise parameters
    m_Tightness = 5;
    m_Alpha = 0.9;
    m_Tempo = 120;
    m_EstimatedTempo = 120.0;
    m_TempoToLagFactor = 60.*44100./512.;
	
    m_m0 = 10;
    m_BeatCounter = -1;
	
    m_BeatDueInFrame = false;
	

	// create rayleigh weighting vector
	for (int n = 0; n < 128; n++)
	{
        m_WeightingVector[n] = ((double) n / pow(rayparam,2)) * exp((-1*pow((double)-n,2)) / (2*pow(rayparam,2)));
	}
	
	// initialise prev_delta
	for (int i = 0; i < 41; i++)
	{
        m_PrevDelta[i] = 1;
	}
	
	double t_mu = 41/2;
	double m_sig;
	double x;
	// create tempo transition matrix
	m_sig = 41/8;
	for (int i = 0;i < 41;i++)
	{
		for (int j = 0;j < 41;j++)
		{
			x = j+1;
			t_mu = i+1;
            m_TempoTransitionMatrix[i][j] = (1 / (m_sig * sqrt(2*pi))) * exp( (-1*pow((x-t_mu),2)) / (2*pow(m_sig,2)) );
		}
	}
	
	// tempo is not fixed
    m_TempoFixed = false;
    
    // initialise latest cumulative score value
    // in case it is requested before any processing takes place
    m_LatestCumulativeScoreValue = 0;
    
    // initialise algorithm given the hopsize and framesize
    updateHopAndFrameSize(hopSize_p, frameSize_p);
    
    
    // Set up FFT for calculating the auto-correlation function
    m_FFTLengthForACFCalculation = 1024;
    
#ifdef USE_FFTW
    m_ComplexIn = (fftw_complex*) fftw_malloc (sizeof(fftw_complex) * m_FFTLengthForACFCalculation);		// complex array to hold fft data
    m_ComplexOut = (fftw_complex*) fftw_malloc (sizeof(fftw_complex) * m_FFTLengthForACFCalculation);	// complex array to hold fft data
    
    //acfForwardFFT = fftw_plan_dft_1d (FFTLengthForACFCalculation, complexIn, complexOut, FFTW_FORWARD, FFTW_ESTIMATE);	// FFT plan initialisation
    //acfBackwardFFT = fftw_plan_dft_1d (FFTLengthForACFCalculation, complexOut, complexIn, FFTW_BACKWARD, FFTW_ESTIMATE);	// FFT plan initialisation
#endif
    
#ifdef USE_KISS_FFT
    m_FFTIn = new kiss_fft_cpx[m_FFTLengthForACFCalculation];
    m_FFTOut = new kiss_fft_cpx[m_FFTLengthForACFCalculation];
    m_CfgForwards = kiss_fft_alloc (m_FFTLengthForACFCalculation, 0, 0, 0);
    m_CfgBackwards = kiss_fft_alloc (m_FFTLengthForACFCalculation, 1, 0, 0);
#endif
}

//=======================================================================
void BTrack::setHopSize (int hopSize_p)
{	
    m_HopSize = hopSize_p;
    m_OnsetDFBufferSize = CIRCULARBUFFERFAC*(512*512)/m_HopSize;		// calculate df buffer size
	
    m_BeatPeriod = round(60/((((double) m_HopSize)/44100)*m_Tempo));

    // set size of onset detection function buffer
    m_OnsetDF.resize (m_OnsetDFBufferSize);
    
    // set size of cumulative score buffer
    m_CumulativeScore.resize (m_OnsetDFBufferSize);

    // set size of magnitude buffer
    m_SpectralPower.resize(m_OnsetDFBufferSize);
	
	// initialise df_buffer to zeros
    for (int i = 0; i < m_OnsetDFBufferSize; i++)
	{
        m_OnsetDF[i] = 0;
        m_CumulativeScore[i] = 0;
		
        if ((i %  ((int) round(m_BeatPeriod))) == 0)
		{
            m_OnsetDF[i] = 1;
		}
	}
}

//=======================================================================
void BTrack::updateHopAndFrameSize (int hopSize_p, int frameSize_p)
{
    // update the onset detection function object
    m_Odf.initialise (hopSize_p, frameSize_p);
    
    // update the hop size being used by the beat tracker
    setHopSize (hopSize_p);
}

//=======================================================================
bool BTrack::beatDueInCurrentFrame()
{
    return m_BeatDueInFrame;
}

//=======================================================================
double BTrack::getCurrentTempoEstimate()
{
    return m_EstimatedTempo;
}

//=======================================================================
int BTrack::getHopSize()
{
    return m_HopSize;
}

//=======================================================================
double BTrack::getLatestCumulativeScoreValue()
{
    return m_LatestCumulativeScoreValue;
}

//=======================================================================
void BTrack::processAudioFrame (double* frame_p)
{
    // calculate the onset detection function sample for the frame
    double sample = m_Odf.calculateOnsetDetectionFunctionSample (frame_p);
    
    // process the new onset detection function sample in the beat tracking algorithm
    processOnsetDetectionFunctionSample (sample);
}

static bool IsStrongPowerTrend(const double timeStep_p, const CircularBuffer& power_p)
{
    using namespace CircularBufferUtils;
    const std::size_t bufferSize = power_p.getCurrentSize();

    if(bufferSize == power_p.getSize())
    {
        /*
        // test
        CircularBuffer testPower;
        testPower.resize(3);
        testPower.addSampleToEnd(0.);
        testPower.addSampleToEnd(1.);
        testPower.addSampleToEnd(2.);
        const LinearFit test = GetLinearFit(timeStep_p, testPower);

        std::cout << test.slope << "\n";*/

        const LinearFit fitResult = GetLinearFit(timeStep_p, power_p);
        const double timeIntervalSize = bufferSize*timeStep_p;
        const double averagePower = std::accumulate(power_p.getBuffer().begin(), power_p.getBuffer().end(), 0.) / bufferSize;

        const double linearChangeEstimate = std::abs(timeIntervalSize*fitResult.slope);
        if(linearChangeEstimate > 10*averagePower)
        {
            return true;
        }
        return false;
    }

    return true; // always trend for not completely filled buffers
}

static bool IsStrongPowerTrend(const CircularBuffer& power_p)
{
    return IsStrongPowerTrend(1.0, power_p); // slope is in units of time step, but time step is multiplied to slope for criterion, hence every time step > 0 is fine for crit
}

//=======================================================================
void BTrack::processOnsetDetectionFunctionSample (double newSample_p)
{
    // we need to ensure that the onset
    // detection function sample is positive
    newSample_p = fabs (newSample_p);
    
    // add a tiny constant to the sample to stop it from ever going
    // to zero. this is to avoid problems further down the line
    newSample_p = newSample_p + 0.0001;
    
    m_m0=std::max(m_m0-1,0);
    m_BeatCounter=std::max(m_BeatCounter-1,0);
    m_BeatDueInFrame = false;
		
	// add new sample at the end
    m_OnsetDF.addSampleToEnd (newSample_p);

    const double magSpecSum = m_Odf.GetMagSpecSum();
    m_SpectralPower.addSampleToEnd(magSpecSum);
	
	// update cumulative score
    updateCumulativeScore (newSample_p);
	
	// if we are halfway between beats
    //bool predictedBeat = false;
    if (m_m0 == 0)
	{
		predictBeat();
        //predictedBeat = true;
	}
	
	// if we are at a beat
    if (m_BeatCounter == 0)
	{
        if(std::abs(newSample_p) > 0.001
                && std::abs(magSpecSum) > 0.3*m_PowerAtLastBeat
                && !IsStrongPowerTrend(m_SpectralPower))
        {
            m_BeatDueInFrame = true;	// indicate a beat should be output
            m_PowerAtLastBeat = std::abs(magSpecSum);
        }

		// recalculate the tempo
		resampleOnsetDetectionFunction();
        calculateTempo(magSpecSum);
	}
}

//=======================================================================
void BTrack::setTempo (double tempo_p)
{	 
	
	/////////// TEMPO INDICATION RESET //////////////////
	

    ReduceTempo(tempo_p);

    m_Tempo = tempo_p;
		
	// convert tempo from bpm value to integer index of tempo probability 
    int tempo_index = (int) round((tempo_p - 80)/2);
	
	// now set previous tempo observations to zero
	for (int i=0;i < 41;i++)
	{
        m_PrevDelta[i] = 0;
	}
	
	// set desired tempo index to 1
    m_PrevDelta[tempo_index] = 1;
	
	
	/////////// CUMULATIVE SCORE ARTIFICAL TEMPO UPDATE //////////////////
	
	// calculate new beat period
    int new_bperiod = (int) round(60/((((double) m_HopSize)/44100)*tempo_p));
	
	int bcounter = 1;
	// initialise df_buffer to zeros
    for (int i = (m_OnsetDFBufferSize-1);i >= 0;i--)
	{
		if (bcounter == 1)
		{
            m_CumulativeScore[i] = 150;
            m_OnsetDF[i] = 150;
		}
		else
		{
            m_CumulativeScore[i] = 10;
            m_OnsetDF[i] = 10;
		}
		
		bcounter++;
		
		if (bcounter > new_bperiod)
		{
			bcounter = 1;
		}
	}
	
	/////////// INDICATE THAT THIS IS A BEAT //////////////////
	
	// beat is now
    m_BeatCounter = 0;
	
	// offbeat is half of new beat period away
    m_m0 = (int) round(((double) new_bperiod)/2);
}

//=======================================================================
void BTrack::fixTempo (double tempo_p)
{	
    // firstly make sure tempo is between 40 and 160 bpm..
    ReduceTempo(tempo_p);

    m_Tempo = tempo_p;
	
	// convert tempo from bpm value to integer index of tempo probability 
    int tempo_index = (int) round((tempo_p - 80)/2);
	
	// now set previous fixed previous tempo observation values to zero
	for (int i=0;i < 41;i++)
	{
        m_PrevDeltaFixed[i] = 0;
	}
	
	// set desired tempo index to 1
    m_PrevDeltaFixed[tempo_index] = 1;
		
	// set the tempo fix flag
    m_TempoFixed = true;
}

//=======================================================================
void BTrack::doNotFixTempo()
{	
	// set the tempo fix flag
    m_TempoFixed = false;
}

//=======================================================================
void BTrack::resampleOnsetDetectionFunction()
{
	float output[512];
    
    std::vector<float>  input(m_OnsetDFBufferSize);
    
    for (int i = 0;i < m_OnsetDFBufferSize;i++)
    {
        input[i] = (float) m_OnsetDF[i];
    }
        
    double src_ratio = 512.0/((double) m_OnsetDFBufferSize);
    int BUFFER_LEN = m_OnsetDFBufferSize;
    int output_len;
    SRC_DATA	src_data ;
    
    //output_len = (int) floor (((double) BUFFER_LEN) * src_ratio) ;
    output_len = 512;
    
    src_data.data_in = input.data();
    src_data.input_frames = BUFFER_LEN;
    
    src_data.src_ratio = src_ratio;
    
    src_data.data_out = output;
    src_data.output_frames = output_len;
    
    src_simple (&src_data, SRC_SINC_BEST_QUALITY, 1);
            
    for (int i = 0;i < output_len;i++)
    {
        m_ResampledOnsetDF[i] = (double) src_data.data_out[i];
    }
}

//=======================================================================
void BTrack::calculateTempo(const double magSpecSum_p)
{
	// adaptive threshold on input
    adaptiveThreshold (m_ResampledOnsetDF,512);
		
	// calculate auto-correlation function of detection function
    calculateBalancedACF (m_ResampledOnsetDF);
	
	// calculate output of comb filterbank
	calculateOutputOfCombFilterBank();
	
	// adaptive threshold on rcf
    adaptiveThreshold (m_CombFilterBankOutput,128);

	
	int t_index;
	int t_index2;
	// calculate tempo observation vector from beat period observation vector
	for (int i = 0;i < 41;i++)
	{
        t_index = (int) round (m_TempoToLagFactor / ((double) ((2*i)+80)));
        t_index2 = (int) round (m_TempoToLagFactor / ((double) ((4*i)+160)));

		
        m_TempoObservationVector[i] = m_CombFilterBankOutput[t_index-1] + m_CombFilterBankOutput[t_index2-1];
	}
	
	
	double maxval;
	double maxind;
	double curval;
	
	// if tempo is fixed then always use a fixed set of tempi as the previous observation probability function
    if (m_TempoFixed)
	{
		for (int k = 0;k < 41;k++)
		{
            m_PrevDelta[k] = m_PrevDeltaFixed[k];
		}
	}
		
	for (int j=0;j < 41;j++)
	{
		maxval = -1;
		for (int i = 0;i < 41;i++)
		{
            curval = m_PrevDelta[i] * m_TempoTransitionMatrix[i][j];
			
			if (curval > maxval)
			{
				maxval = curval;
			}
		}
		
        m_Delta[j] = maxval * m_TempoObservationVector[j];
	}
	

    normaliseArray(m_Delta,41);
	
	maxind = -1;
	maxval = -1;
	
	for (int j=0;j < 41;j++)
	{
        if (m_Delta[j] > maxval)
		{
            maxval = m_Delta[j];
			maxind = j;
		}
		
        m_PrevDelta[j] = m_Delta[j];
	}
	
    m_BeatPeriod = round ((60.0*44100.0)/(((2*maxind)+80)*((double) m_HopSize)));
	
    if (m_BeatPeriod > 0)
	{
        if(std::abs(magSpecSum_p) < 0.3*m_PowerAtLastBeat || std::abs(magSpecSum_p) < 0.1)
        {
            m_EstimatedTempo *= 0.9;
        }
        else
        {
            m_EstimatedTempo = 60.0/((((double) m_HopSize) / 44100.0) * m_BeatPeriod);
        }

        ReduceTempo(m_EstimatedTempo);
	}
}

//=======================================================================
void BTrack::adaptiveThreshold (double* x_p, int N_p)
{
	int i = 0;
	int k,t = 0;
    std::vector<double>  x_thresh(N_p);
	
	int p_post = 7;
	int p_pre = 8;
	
    t = std::min(N_p,p_post);	// what is smaller, p_post of df size. This is to avoid accessing outside of arrays
	
	// find threshold for first 't' samples, where a full average cannot be computed yet 
	for (i = 0;i <= t;i++)
	{	
        k = std::min ((i+p_pre),N_p);
        x_thresh[i] = calculateMeanOfArray (x_p,1,k);
	}
	// find threshold for bulk of samples across a moving average from [i-p_pre,i+p_post]
    for (i = t+1;i < N_p-p_post;i++)
	{
        x_thresh[i] = calculateMeanOfArray (x_p,i-p_pre,i+p_post);
	}
	// for last few samples calculate threshold, again, not enough samples to do as above
    for (i = N_p-p_post;i < N_p;i++)
	{
		k = std::max ((i-p_post),1);
        x_thresh[i] = calculateMeanOfArray (x_p,k,N_p);
	}
	
	// subtract the threshold from the detection function and check that it is not less than 0
    for (i = 0; i < N_p; i++)
	{
        x_p[i] = x_p[i] - x_thresh[i];
        if (x_p[i] < 0)
		{
            x_p[i] = 0;
		}
	}
}

//=======================================================================
void BTrack::calculateOutputOfCombFilterBank()
{
	int numelem;
	
	for (int i = 0;i < 128;i++)
	{
        m_CombFilterBankOutput[i] = 0;
	}
	
	numelem = 4;
	
	for (int i = 2; i <= 127; i++) // max beat period
	{
		for (int a = 1; a <= numelem; a++) // number of comb elements
		{
			for (int b = 1-a; b <= a-1; b++) // general state using normalisation of comb elements
			{
                m_CombFilterBankOutput[i-1] = m_CombFilterBankOutput[i-1] + (m_Acf[(a*i+b)-1]*m_WeightingVector[i-1])/(2*a-1);	// calculate value for comb filter row
			}
		}
	}
}

//=======================================================================
void BTrack::calculateBalancedACF (double* onsetDetectionFunction_p)
{
    int onsetDetectionFunctionLength = 512;
    
#ifdef USE_FFTW
    // copy into complex array and zero pad
    for (int i = 0;i < m_FFTLengthForACFCalculation;i++)
    {
        if (i < onsetDetectionFunctionLength)
        {
            m_ComplexIn[i][0] = onsetDetectionFunction_p[i];
            m_ComplexIn[i][1] = 0.0;
        }
        else
        {
            m_ComplexIn[i][0] = 0.0;
            m_ComplexIn[i][1] = 0.0;
        }
    }
    
    // perform the fft
    fftw_plan acfForwardFFT = fftw_plan_dft_1d (m_FFTLengthForACFCalculation, m_ComplexIn, m_ComplexOut, FFTW_FORWARD, 0);
    fftw_execute(acfForwardFFT);
    fftw_destroy_plan(acfForwardFFT);
    
    // multiply by complex conjugate
    for (int i = 0;i < m_FFTLengthForACFCalculation;i++)
    {
        m_ComplexOut[i][0] = m_ComplexOut[i][0]*m_ComplexOut[i][0] + m_ComplexOut[i][1]*m_ComplexOut[i][1];
        m_ComplexOut[i][1] = 0.0;
    }
    
    // perform the ifft
    fftw_plan acfBackwardFFT = fftw_plan_dft_1d (m_FFTLengthForACFCalculation, m_ComplexOut, m_ComplexIn, FFTW_BACKWARD, 0);
    fftw_execute (acfBackwardFFT);
    fftw_destroy_plan(acfBackwardFFT);
    
#endif
    
#ifdef USE_KISS_FFT
    // copy into complex array and zero pad
    for (int i = 0;i < m_FFTLengthForACFCalculation;i++)
    {
        if (i < onsetDetectionFunctionLength)
        {
            m_FFTIn[i].r = onsetDetectionFunction_p[i];
            m_FFTIn[i].i = 0.0;
        }
        else
        {
            m_FFTIn[i].r = 0.0;
            m_FFTIn[i].i = 0.0;
        }
    }
    
    // execute kiss fft
    kiss_fft (m_CfgForwards, m_FFTIn, m_FFTOut);
    
    // multiply by complex conjugate
    for (int i = 0;i < m_FFTLengthForACFCalculation;i++)
    {
        m_FFTOut[i].r = m_FFTOut[i].r * m_FFTOut[i].r + m_FFTOut[i].i * m_FFTOut[i].i;
        m_FFTOut[i].i = 0.0;
    }
    
    // perform the ifft
    kiss_fft (m_CfgBackwards, m_FFTOut, m_FFTIn);
    
#endif
    
    double lag = 512;
    
    for (int i = 0; i < 512; i++)
    {
#ifdef USE_FFTW
        // calculate absolute value of result
        double absValue = sqrt (m_ComplexIn[i][0]*m_ComplexIn[i][0] + m_ComplexIn[i][1]*m_ComplexIn[i][1]);
#endif
        
#ifdef USE_KISS_FFT
        // calculate absolute value of result
        double absValue = sqrt (m_FFTIn[i].r * m_FFTIn[i].r + m_FFTIn[i].i * m_FFTIn[i].i);
#endif
        // divide by inverse lad to deal with scale bias towards small lags
        m_Acf[i] = absValue / lag;
        
        // this division by 1024 is technically unnecessary but it ensures the algorithm produces
        // exactly the same ACF output as the old time domain implementation. The time difference is
        // minimal so I decided to keep it
        m_Acf[i] = m_Acf[i] / 1024.;
        
        lag = lag - 1.;
    }
}

//=======================================================================
double BTrack::calculateMeanOfArray (double* array_p, int startIndex_p, int endIndex_p)
{
	int i;
	double sum = 0;

    int length = endIndex_p - startIndex_p;
	
	// find sum
    for (i = startIndex_p; i < endIndex_p; i++)
	{
        sum = sum + array_p[i];
	}
	
    if (length > 0)
    {
        return sum / length;	// average and return
    }
    else
    {
        return 0;
    }
}

//=======================================================================
void BTrack::normaliseArray (double* array_p, int N_p)
{
	double sum = 0;
	
    for (int i = 0; i < N_p; i++)
	{
        if (array_p[i] > 0)
		{
            sum = sum + array_p[i];
		}
	}
	
	if (sum > 0)
	{
        for (int i = 0; i < N_p; i++)
		{
            array_p[i] = array_p[i] / sum;
		}
	}
}

//=======================================================================
void BTrack::updateCumulativeScore (double odfSample_p)
{
	int start, end, winsize;
	double max;
	
    start = m_OnsetDFBufferSize - round (2 * m_BeatPeriod);
    end = m_OnsetDFBufferSize - round (m_BeatPeriod / 2);
	winsize = end-start+1;
	
    std::vector<double> w1(winsize);
    double v = -2*m_BeatPeriod;
	double wcumscore;
	
	// create window
	for (int i = 0; i < winsize; i++)
	{
        w1[i] = exp((-1 * pow (m_Tightness * log (-v / m_BeatPeriod), 2)) / 2);
		v = v+1;
	}	
	
	// calculate new cumulative score value
	max = 0;
	int n = 0;
    for (int i=start; i <= end; i++)
    {
        wcumscore = m_CumulativeScore[i]*w1[n];

        if (wcumscore > max)
        {
            max = wcumscore;
        }
        n++;
    }
	
    m_LatestCumulativeScoreValue = ((1 - m_Alpha) * odfSample_p) + (m_Alpha * max);
    
    m_CumulativeScore.addSampleToEnd (m_LatestCumulativeScoreValue);
}

//=======================================================================
void BTrack::predictBeat()
{	 
    int windowSize = (int) m_BeatPeriod;
    std::vector<double> futureCumulativeScore(m_OnsetDFBufferSize + windowSize);
    std::vector<double> w2(windowSize);
    
	// copy cumscore to first part of fcumscore
    for (int i = 0;i < m_OnsetDFBufferSize;i++)
	{
        futureCumulativeScore[i] = m_CumulativeScore[i];
	}
	
	// create future window
	double v = 1;
	for (int i = 0; i < windowSize; i++)
	{
        w2[i] = exp((-1*pow((v - (m_BeatPeriod/2)),2))   /  (2*pow((m_BeatPeriod/2) ,2)));
		v++;
	}
	
	// create past window
    v = -2*m_BeatPeriod;
    int start = m_OnsetDFBufferSize - round(2*m_BeatPeriod);
    int end = m_OnsetDFBufferSize - round(m_BeatPeriod/2);
	int pastwinsize = end-start+1;
    std::vector<double> w1(pastwinsize);

	for (int i = 0;i < pastwinsize;i++)
	{
        w1[i] = exp((-1*pow(m_Tightness*log(-v/m_BeatPeriod),2))/2);
		v = v+1;
	}

	// calculate future cumulative score
	double max;
	int n;
	double wcumscore;
    for (int i = m_OnsetDFBufferSize; i < (m_OnsetDFBufferSize + windowSize); i++)
	{
        start = i - round (2*m_BeatPeriod);
        end = i - round (m_BeatPeriod/2);
		
		max = 0;
		n = 0;
		for (int k=start;k <= end;k++)
		{
			wcumscore = futureCumulativeScore[k]*w1[n];
			
			if (wcumscore > max)
			{
				max = wcumscore;
			}
			n++;
		}
		
		futureCumulativeScore[i] = max;
	}
	
	// predict beat
	max = 0;
	n = 0;
	
    for (int i = m_OnsetDFBufferSize; i < (m_OnsetDFBufferSize + windowSize); i++)
	{
		wcumscore = futureCumulativeScore[i]*w2[n];
		
		if (wcumscore > max)
		{
			max = wcumscore;
            m_BeatCounter = n;
		}	
		
		n++;
	}
		
	// set next prediction time
    m_m0 = m_BeatCounter + round (m_BeatPeriod / 2);
}
