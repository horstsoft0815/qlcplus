//=======================================================================
/** @file OnsetDetectionFunction.cpp
 *  @brief A class for calculating onset detection functions
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

#include <math.h>
#include "OnsetDetectionFunction.h"

#define CLEAR_FFT_NOISE

//=======================================================================
OnsetDetectionFunction::OnsetDetectionFunction (
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
OnsetDetectionFunction::OnsetDetectionFunction(
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
OnsetDetectionFunction::~OnsetDetectionFunction()
{
    if (m_Initialised)
    {
        //freeFFT();
    }
}

//=======================================================================
void OnsetDetectionFunction::initialise (int hopSize_p, int frameSize_p)
{
    // use the already initialised onset detection function and window type and
    // pass the new frame and hop size to the main initialisation function
    initialise (hopSize_p, frameSize_p, m_OnsetDetectionFunctionType, m_WindowType);
}

//=======================================================================
void OnsetDetectionFunction::initialise(
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
void OnsetDetectionFunction::initialiseFFT()
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
    

    m_FFTIn = new kiss_fft_cpx[m_FrameSize];
    m_FFTOut = new kiss_fft_cpx[m_FrameSize];
    m_Cfg = kiss_fft_alloc (m_FrameSize, 0, 0, 0);
#endif

    m_Initialised = true;
}

//=======================================================================
void OnsetDetectionFunction::freeFFT()
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

//=======================================================================
void OnsetDetectionFunction::setOnsetDetectionFunctionType (
        int onsetDetectionFunctionType_p)
{
    m_OnsetDetectionFunctionType = onsetDetectionFunctionType_p; // set detection function type
}

//=======================================================================
double OnsetDetectionFunction::calculateOnsetDetectionFunctionSample (double* buffer_p)
{	
	double odfSample;
		
	// shift audio samples back in frame by hop size
    for (int i = 0; i < (m_FrameSize-m_HopSize);i++)
	{
        m_Frame[i] = m_Frame[i+m_HopSize];
	}
	
	// add new samples to frame from input buffer
	int j = 0;
    for (int i = (m_FrameSize-m_HopSize);i < m_FrameSize;i++)
	{
        m_Frame[i] = buffer_p[j];
		j++;
	}
		
    switch (m_OnsetDetectionFunctionType)
    {
        case eEnergyEnvelope:
        {
            // calculate energy envelope detection function sample
			odfSample = energyEnvelope();
			break;
        }
        case eEnergyDifference:
        {
            // calculate half-wave rectified energy difference detection function sample
			odfSample = energyDifference();
			break;
        }
        case eSpectralDifference:
        {
            // calculate spectral difference detection function sample
			odfSample = spectralDifference();
			break;
        }
        case eSpectralDifferenceHWR:
        {
            // calculate spectral difference detection function sample (half wave rectified)
			odfSample = spectralDifferenceHWR();
			break;
        }
        case ePhaseDeviation:
        {
            // calculate phase deviation detection function sample (half wave rectified)
			odfSample = phaseDeviation();
			break;
        }
        case eComplexSpectralDifference:
        {
            // calcualte complex spectral difference detection function sample
			odfSample = complexSpectralDifference();
			break;
        }
        case eComplexSpectralDifferenceHWR:
        {
            // calcualte complex spectral difference detection function sample (half-wave rectified)
			odfSample = complexSpectralDifferenceHWR();
			break;
        }
        case eHighFrequencyContent:
        {
            // calculate high frequency content detection function sample
			odfSample = highFrequencyContent();
			break;
        }
        case eHighFrequencySpectralDifference:
        {
            // calculate high frequency spectral difference detection function sample
			odfSample = highFrequencySpectralDifference();
			break;
        }
        case eHighFrequencySpectralDifferenceHWR:
        {
            // calculate high frequency spectral difference detection function (half-wave rectified)
			odfSample = highFrequencySpectralDifferenceHWR();
			break;
        }
		default:
        {
			odfSample = 1.0;
        }
	}
		
    return odfSample;
}


//=======================================================================
void OnsetDetectionFunction::performFFT()
{
    int fsize2 = (m_FrameSize/2);
    
#ifdef USE_FFTW
	// window frame and copy to complex array, swapping the first and second half of the signal
	for (int i = 0;i < fsize2;i++)
	{
        m_ComplexIn[i][0] = m_Frame[i + fsize2] * m_Window[i + fsize2];
        m_ComplexIn[i][1] = 0.0;
        m_ComplexIn[i+fsize2][0] = m_Frame[i] * m_Window[i];
        m_ComplexIn[i+fsize2][1] = 0.0;
	}
	
	// perform the fft
    fftw_plan m_Plan = fftw_plan_dft_1d (m_FrameSize, m_ComplexIn, m_ComplexOut, FFTW_FORWARD, FFTW_ESTIMATE);	// FFT plan initialisation
    fftw_execute (m_Plan);
    fftw_destroy_plan (m_Plan);
#endif
    
#ifdef USE_KISS_FFT
    for (int i = 0; i < fsize2; i++)
    {
        m_FFTIn[i].r = m_Frame[i + fsize2] * m_Window[i + fsize2];
        m_FFTIn[i].i = 0.0;
        m_FFTIn[i + fsize2].r = m_Frame[i] * m_Window[i];
        m_FFTIn[i + fsize2].i = 0.0;
    }
    
    // execute kiss fft
    kiss_fft (m_Cfg, m_FFTIn, m_FFTOut);
    
    // store real and imaginary parts of FFT
    for (int i = 0; i < m_FrameSize; i++)
    {
        m_ComplexOut[i][0] = m_FFTOut[i].r;
        m_ComplexOut[i][1] = m_FFTOut[i].i;
    }
#endif

#ifdef CLEAR_FFT_NOISE
    //We delete some values since these will ruin our output
    for (size_t n = 0; n < 5 && n < m_ComplexOut.size(); n++)
    {
        m_ComplexOut[n][0] = 0;
        m_ComplexOut[n][1] = 0;
    }

    // consider only data < 5000Hz
    double SPECTRUM_MAX_FREQUENCY = 2500;
    double sampleRate = 44100;
    int number_of_bands = 1;
    size_t subBandWidth = static_cast<size_t>(((m_ComplexOut.size() * SPECTRUM_MAX_FREQUENCY) / sampleRate) / number_of_bands);

    for (size_t n = subBandWidth; n < m_ComplexOut.size(); n++)
    {
        m_ComplexOut[n][0] = 0;
        m_ComplexOut[n][1] = 0;
    }
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Methods for Detection Functions /////////////////////////////////

//=======================================================================
double OnsetDetectionFunction::energyEnvelope()
{
	double sum;
	
	sum = 0;	// initialise sum
	
	// sum the squares of the samples
    for (int i = 0;i < m_FrameSize;i++)
	{
        sum = sum + (m_Frame[i] * m_Frame[i]);
	}
	
	return sum;		// return sum
}

//=======================================================================
double OnsetDetectionFunction::energyDifference()
{
	double sum;
	double sample;
	
	sum = 0;	// initialise sum
	
	// sum the squares of the samples
    for (int i = 0; i < m_FrameSize; i++)
	{
        sum = sum + (m_Frame[i] * m_Frame[i]);
	}
	
    sample = sum - m_PrevEnergySum;	// sample is first order difference in energy
	
    m_PrevEnergySum = sum;	// store energy value for next calculation
	
	if (sample > 0)
	{
		return sample;		// return difference
	}
	else
	{
		return 0;
	}
}

//=======================================================================
double OnsetDetectionFunction::spectralDifference()
{
	double diff;
	double sum;
	
	// perform the FFT
	performFFT();
	
	// compute first (N/2)+1 mag values
    for (int i = 0;i < (m_FrameSize/2)+1;i++)
	{
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0], 2) + pow (m_ComplexOut[i][1], 2));
	}
	// mag spec symmetric above (N/2)+1 so copy previous values
    for (int i = (m_FrameSize/2)+1; i < m_FrameSize; i++)
	{
        m_MagSpec[i] = m_MagSpec[m_FrameSize-i];
	}
	
	sum = 0;	// initialise sum to zero
    m_MagSpecSum = 0;
    for (int i = 0; i < m_FrameSize; i++)
	{
		// calculate difference
        diff = m_MagSpec[i] - m_PrevMagSpec[i];
        m_MagSpecSum += m_MagSpec[i];

		// ensure all difference values are positive
		if (diff < 0)
		{
			diff = diff*-1;
		}
		
		// add difference to sum
		sum = sum + diff;
		
		// store magnitude spectrum bin for next detection function sample calculation
        m_PrevMagSpec[i] = m_MagSpec[i];
	}
	
	return sum;		
}

//=======================================================================
double OnsetDetectionFunction::spectralDifferenceHWR()
{
	double diff;
	double sum;
	
	// perform the FFT
	performFFT();
	
	// compute first (N/2)+1 mag values
    for (int i = 0;i < (m_FrameSize/2) + 1; i++)
	{
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0],2) + pow (m_ComplexOut[i][1],2));
	}
	// mag spec symmetric above (N/2)+1 so copy previous values
    for (int i = (m_FrameSize/2)+1;i < m_FrameSize;i++)
	{
        m_MagSpec[i] = m_MagSpec[m_FrameSize-i];
	}
	
	sum = 0;	// initialise sum to zero
    m_MagSpecSum = 0;
    for (int i = 0;i < m_FrameSize;i++)
	{
		// calculate difference
        diff = m_MagSpec[i] - m_PrevMagSpec[i];
        m_MagSpecSum += m_MagSpec[i];

		// only add up positive differences
		if (diff > 0)
		{
			// add difference to sum
			sum = sum+diff;
		}
		
		// store magnitude spectrum bin for next detection function sample calculation
        m_PrevMagSpec[i] = m_MagSpec[i];
	}
	
	return sum;		
}


//=======================================================================
double OnsetDetectionFunction::phaseDeviation()
{
	double dev,pdev;
	double sum;
	
	// perform the FFT
	performFFT();
	
	sum = 0; // initialise sum to zero
    m_MagSpecSum = 0;
	// compute phase values from fft output and sum deviations
    for (int i = 0;i < m_FrameSize;i++)
	{
		// calculate phase value
        m_Phase[i] = atan2 (m_ComplexOut[i][1], m_ComplexOut[i][0]);
		
		// calculate magnitude value
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0],2) + pow (m_ComplexOut[i][1],2));
        m_MagSpecSum += m_MagSpec[i];
		
		// if bin is not just a low energy bin then examine phase deviation
        if (m_MagSpec[i] > 0.1)
		{
            dev = m_Phase[i] - (2*m_PrevPhase[i]) + m_PrevPhase2[i];	// phase deviation
			pdev = princarg (dev);	// wrap into [-pi,pi] range
		
			// make all values positive
			if (pdev < 0)	
			{
				pdev = pdev*-1;
			}
						
			// add to sum
			sum = sum + pdev;
		}
				
		// store values for next calculation
        m_PrevPhase2[i] = m_PrevPhase[i];
        m_PrevPhase[i] = m_Phase[i];
	}
	
	return sum;		
}

//=======================================================================
double OnsetDetectionFunction::complexSpectralDifference()
{
	double phaseDeviation;
	double sum;
	double csd;
	
	// perform the FFT
	performFFT();
	
	sum = 0; // initialise sum to zero
    m_MagSpecSum = 0;
	// compute phase values from fft output and sum deviations
    for (int i = 0;i < m_FrameSize;i++)
	{
		// calculate phase value
        m_Phase[i] = atan2 (m_ComplexOut[i][1], m_ComplexOut[i][0]);
		
		// calculate magnitude value
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0],2) + pow(m_ComplexOut[i][1],2));
        m_MagSpecSum += m_MagSpec[i];
		
		// phase deviation
        phaseDeviation = m_Phase[i] - (2 * m_PrevPhase[i]) + m_PrevPhase2[i];
		
        // calculate complex spectral difference for the current spectral bin
        csd = sqrt (pow (m_MagSpec[i], 2) + pow (m_PrevMagSpec[i], 2) - 2 * m_MagSpec[i] * m_PrevMagSpec[i] * cos (phaseDeviation));
			
		// add to sum
		sum = sum + csd;
		
		// store values for next calculation
        m_PrevPhase2[i] = m_PrevPhase[i];
        m_PrevPhase[i] = m_Phase[i];
        m_PrevMagSpec[i] = m_MagSpec[i];
	}
	
	return sum;		
}

//=======================================================================
double OnsetDetectionFunction::complexSpectralDifferenceHWR()
{
	double phaseDeviation;
	double sum;
	double magnitudeDifference;
	double csd;
	
	// perform the FFT
	performFFT();
	
	sum = 0; // initialise sum to zero
    m_MagSpecSum = 0;
	// compute phase values from fft output and sum deviations
    for (int i = 0;i < m_FrameSize;i++)
	{
		// calculate phase value
        m_Phase[i] = atan2 (m_ComplexOut[i][1], m_ComplexOut[i][0]);
		
		// calculate magnitude value
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0],2) + pow(m_ComplexOut[i][1],2));
        m_MagSpecSum += m_MagSpec[i];
		
        // phase deviation
        phaseDeviation = m_Phase[i] - (2 * m_PrevPhase[i]) + m_PrevPhase2[i];
        
        // calculate magnitude difference (real part of Euclidean distance between complex frames)
        magnitudeDifference = m_MagSpec[i] - m_PrevMagSpec[i];
        
        // if we have a positive change in magnitude, then include in sum, otherwise ignore (half-wave rectification)
        if (magnitudeDifference > 0)
        {
            // calculate complex spectral difference for the current spectral bin
            csd = sqrt (pow (m_MagSpec[i], 2) + pow (m_PrevMagSpec[i], 2) - 2 * m_MagSpec[i] * m_PrevMagSpec[i] * cos (phaseDeviation));
        
            // add to sum
            sum = sum + csd;
        }
        
		// store values for next calculation
        m_PrevPhase2[i] = m_PrevPhase[i];
        m_PrevPhase[i] = m_Phase[i];
        m_PrevMagSpec[i] = m_MagSpec[i];
	}
	
	return sum;		
}


//=======================================================================
double OnsetDetectionFunction::highFrequencyContent()
{
	double sum;
	
	// perform the FFT
	performFFT();
	
	sum = 0; // initialise sum to zero
    m_MagSpecSum = 0;
	// compute phase values from fft output and sum deviations
    for (int i = 0; i < m_FrameSize; i++)
	{		
		// calculate magnitude value
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0],2) + pow (m_ComplexOut[i][1],2));
        m_MagSpecSum += m_MagSpec[i];
		
        sum = sum + (m_MagSpec[i] * ((double) (i+1)));
		
		// store values for next calculation
        m_PrevMagSpec[i] = m_MagSpec[i];
	}
	
	return sum;		
}

//=======================================================================
double OnsetDetectionFunction::highFrequencySpectralDifference()
{
	double sum;
	double mag_diff;
	
	// perform the FFT
	performFFT();
	
	sum = 0; // initialise sum to zero
    m_MagSpecSum = 0;
	// compute phase values from fft output and sum deviations
    for (int i = 0;i < m_FrameSize;i++)
	{		
		// calculate magnitude value
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0],2) + pow (m_ComplexOut[i][1],2));
        m_MagSpecSum += m_MagSpec[i];

		// calculate difference
        mag_diff = m_MagSpec[i] - m_PrevMagSpec[i];
		
		if (mag_diff < 0)
		{
			mag_diff = -mag_diff;
		}
		
		sum = sum + (mag_diff * ((double) (i+1)));
		
		// store values for next calculation
        m_PrevMagSpec[i] = m_MagSpec[i];
	}
	
	return sum;		
}

//=======================================================================
double OnsetDetectionFunction::highFrequencySpectralDifferenceHWR()
{
	double sum;
	double mag_diff;
	
	// perform the FFT
	performFFT();
	
	sum = 0; // initialise sum to zero
    m_MagSpecSum = 0;
	// compute phase values from fft output and sum deviations
    for (int i = 0;i < m_FrameSize;i++)
	{		
		// calculate magnitude value
        m_MagSpec[i] = sqrt (pow (m_ComplexOut[i][0],2) + pow (m_ComplexOut[i][1],2));
        m_MagSpecSum += m_MagSpec[i];
		
		// calculate difference
        mag_diff = m_MagSpec[i] - m_PrevMagSpec[i];
		
		if (mag_diff > 0)
		{
			sum = sum + (mag_diff * ((double) (i+1)));
		}

		// store values for next calculation
        m_PrevMagSpec[i] = m_MagSpec[i];
	}
	
	return sum;		
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Methods to Calculate Windows ////////////////////////////////////

//=======================================================================
void OnsetDetectionFunction::calculateHanningWindow()
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
void OnsetDetectionFunction::calclulateHammingWindow()
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
void OnsetDetectionFunction::calculateBlackmanWindow()
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
void OnsetDetectionFunction::calculateTukeyWindow()
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
void OnsetDetectionFunction::calculateRectangularWindow()
{
	// Rectangular window calculation
    for (int n = 0;n < m_FrameSize;n++)
	{
        m_Window[n] = 1.0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Other Handy Methods //////////////////////////////////////////

//=======================================================================
double OnsetDetectionFunction::princarg(double phaseVal)
{	
	// if phase value is less than or equal to -pi then add 2*pi
    while (phaseVal <= (-m_Pi))
	{
        phaseVal = phaseVal + (2 * m_Pi);
	}
	
	// if phase value is larger than pi, then subtract 2*pi
    while (phaseVal > m_Pi)
	{
        phaseVal = phaseVal - (2 * m_Pi);
	}
			
	return phaseVal;
}
