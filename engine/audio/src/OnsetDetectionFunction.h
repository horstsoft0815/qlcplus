//=======================================================================
/** @file OnsetDetectionFunction.h
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

#ifndef __ONSETDETECTIONFUNCTION_H
#define __ONSETDETECTIONFUNCTION_H

#ifdef USE_FFTW
#include "fftw3.h"
#endif

#ifdef USE_KISS_FFT
#include "kiss_fft.h"
#endif

#include <vector>
#include <array>
//=======================================================================
/** The type of onset detection function to calculate */
enum OnsetDetectionFunctionType
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

//=======================================================================
/** The type of window to use when calculating onset detection function samples */
enum WindowType
{
    eRectangularWindow,
    eHanningWindow,
    eHammingWindow,
    eBlackmanWindow,
    eTukeyWindow
};

//=======================================================================
/** A class for calculating onset detection functions. */
class OnsetDetectionFunction
{
public:
    
    /** Constructor that defaults the onset detection function type to ComplexSpectralDifferenceHWR
     * and the window type to HanningWindow
     * @param hopSize_p the hop size in audio samples
     * @param frameSize_p the frame size in audio samples
     */
    OnsetDetectionFunction (int hopSize_p, int frameSize_p);
    
    
    /** Constructor 
     * @param hopSize_p the hop size in audio samples
     * @param frameSize_p the frame size in audio samples
     * @param onsetDetectionFunctionType_p the type of onset detection function to use - (see OnsetDetectionFunctionType)
     * @param windowType_p the type of window to use (see WindowType)
     */
    OnsetDetectionFunction (int hopSize_p, int frameSize_p, int onsetDetectionFunctionType_p, int windowType_p);
    
    /** Destructor */
	~OnsetDetectionFunction();
    
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
	
    /** Process input frame and calculate detection function sample 
     * @param buffer_p a pointer to an array containing the audio samples to be processed
     * @returns the onset detection function sample
     */
    double calculateOnsetDetectionFunctionSample (double* buffer_p);
    
    /** Set the detection function type 
     * @param onsetDetectionFunctionType_p the type of onset detection function to use - (see OnsetDetectionFunctionType)
     */
    void setOnsetDetectionFunctionType (int onsetDetectionFunctionType_p);

    double GetMagSpecSum() const { return m_MagSpecSum; }
	
private:
	
    /** Perform the FFT on the data in 'frame' */
	void performFFT();

    //=======================================================================
    /** Calculate energy envelope detection function sample */
	double energyEnvelope();
    
    /** Calculate energy difference detection function sample */
	double energyDifference();
    
    /** Calculate spectral difference detection function sample */
	double spectralDifference();
    
    /** Calculate spectral difference (half wave rectified) detection function sample */
	double spectralDifferenceHWR();
    
    /** Calculate phase deviation detection function sample */
	double phaseDeviation();
    
    /** Calculate complex spectral difference detection function sample */
	double complexSpectralDifference();
    
    /** Calculate complex spectral difference detection function sample (half-wave rectified) */
	double complexSpectralDifferenceHWR();
    
    /** Calculate high frequency content detection function sample */
	double highFrequencyContent();
    
    /** Calculate high frequency spectral difference detection function sample */
	double highFrequencySpectralDifference();
    
    /** Calculate high frequency spectral difference detection function sample (half-wave rectified) */
	double highFrequencySpectralDifferenceHWR();

    //=======================================================================
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

    //=======================================================================
	/** Set phase values between [-pi, pi] 
     * @param phaseVal_p the phase value to process
     * @returns the wrapped phase value
     */
    double princarg(double phaseVal_p);
	
    void initialiseFFT();
    void freeFFT();
	
    double m_Pi;							/**< pi, the constant */
	
    int m_FrameSize;						/**< audio framesize */
    int m_HopSize;						/**< audio hopsize */
    int m_OnsetDetectionFunctionType;		/**< type of detection function */
    int m_WindowType;                     /**< type of window used in calculations */

    //=======================================================================
#ifdef USE_FFTW
    //fftw_plan m_Plan;						/**< fftw plan */
    fftw_complex* m_ComplexIn;			/**< to hold complex fft values for input */
    fftw_complex* m_ComplexOut;			/**< to hold complex fft values for output */
#endif
    
#ifdef USE_KISS_FFT
    kiss_fft_cfg m_Cfg;                   /**< Kiss FFT configuration */
    kiss_fft_cpx* m_FFTIn;                /**< FFT input samples, in complex form */
    kiss_fft_cpx* m_FFTOut;               /**< FFT output samples, in complex form */
    std::vector<std::array<double,2> >  m_ComplexOut;
#endif
	
    //=======================================================================
    bool m_Initialised = false;					/**< flag indicating whether buffers and FFT plans are initialised */

    std::vector<double> m_Frame;          /**< audio frame */
    std::vector<double> m_Window;         /**< window */
	
    double m_PrevEnergySum;				/**< to hold the previous energy sum value */
	
    double m_MagSpecSum = 0;

    std::vector<double> m_MagSpec;        /**< magnitude spectrum */
    std::vector<double> m_PrevMagSpec;    /**< previous magnitude spectrum */
	
    std::vector<double> m_Phase;          /**< FFT phase values */
    std::vector<double> m_PrevPhase;      /**< previous phase values */
    std::vector<double> m_PrevPhase2;     /**< second order previous phase values */

};


#endif
