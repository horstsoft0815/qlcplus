//=======================================================================
/** @file BTrack.h
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

#ifndef __BTRACK_H
#define __BTRACK_H

#include "OnsetDetectionFunction.h"
#include "CircularBuffer.h"
#include <vector>

//=======================================================================
/** The main beat tracking class and the interface to the BTrack
 * beat tracking algorithm. The algorithm can process either
 * audio frames or onset detection function samples and also
 * contains some static functions for calculating beat times in seconds
 */
class BTrack {
	
public:
    
    //=======================================================================
    /** Constructor assuming hop size of 512 and frame size of 1024 */
    BTrack();
    
    /** Constructor assuming frame size will be double the hopSize
     * @param hopSize_p the hop size in audio samples
     */
    BTrack (int hopSize_p);
    
    /** Constructor taking both hopSize and frameSize
     * @param hopSize_p the hop size in audio samples
     * @param frameSize_p the frame size in audio samples
     */
    BTrack (int hopSize_p, int frameSize_p);
    
    /** Destructor */
    ~BTrack();
    
    //=======================================================================
    /** Updates the hop and frame size used by the beat tracker 
     * @param hopSize_p the hop size in audio samples
     * @param frameSize_p the frame size in audio samples
     */
    void updateHopAndFrameSize (int hopSize_p, int frameSize_p);
    
    //=======================================================================
    /** Process a single audio frame 
     * @param frame_p a pointer to an array containing an audio frame. The number of samples should
     * match the frame size that the algorithm was initialised with.
     */
    void processAudioFrame (double* frame_p);
    
    /** Add new onset detection function sample to buffer and apply beat tracking 
     * @param sample_p an onset detection function sample
     */
    void processOnsetDetectionFunctionSample (double sample_p/*, double newMagnitudeSum_p*/);
   
    //=======================================================================
    /** @returns the current hop size being used by the beat tracker */
    int getHopSize();
    
    /** @returns true if a beat should occur in the current audio frame */
    bool beatDueInCurrentFrame();

    /** @returns the current tempo estimate being used by the beat tracker */
    double getCurrentTempoEstimate();
    
    /** @returns the most recent value of the cumulative score function */
    double getLatestCumulativeScoreValue();
    
    //=======================================================================
    /** Set the tempo of the beat tracker 
     * @param tempo_p the tempo in beats per minute (bpm)
     */
    void setTempo (double tempo_p);
    double getTempo () const { return m_Tempo; }
    
    /** Fix tempo to roughly around some value, so that the algorithm will only try to track
     * tempi around the given tempo
     * @param tempo_p the tempo in beats per minute (bpm)
     */
    void fixTempo (double tempo_p);
    
    /** Tell the algorithm to not fix the tempo anymore */
    void doNotFixTempo();
    
    //=======================================================================
    /** Calculates a beat time in seconds, given the frame number, hop size and sampling frequency.
     * This version uses a long to represent the frame number
     * @param frameNumber_p the index of the current frame
     * @param hopSize_p the hop size in audio samples
     * @param fs_p the sampling frequency in Hz
     * @returns a beat time in seconds
     */
    static double getBeatTimeInSeconds (long frameNumber_p, int hopSize_p, int fs_p);
    
    /** Calculates a beat time in seconds, given the frame number, hop size and sampling frequency.
     * This version uses an int to represent the frame number
     * @param frameNumber_p the index of the current frame
     * @param hopSize_p the hop size in audio samples
     * @param fs_p the sampling frequency in Hz
     * @returns a beat time in seconds
     */
    static double getBeatTimeInSeconds (int frameNumber_p, int hopSize_p, int fs_p);
    
		
private:
    
    /** Initialises the algorithm, setting internal parameters and creating weighting vectors 
     * @param hopSize_p the hop size in audio samples
     * @param frameSize_p the frame size in audio samples
     */
    void initialise (int hopSize_p, int frameSize_p);
    
    /** Initialise with hop size and set all array sizes accordingly
     * @param hopSize_p the hop size in audio samples
     */
    void setHopSize (int hopSize_p);
    
    /** Resamples the onset detection function from an arbitrary number of samples to 512 */
    void resampleOnsetDetectionFunction();
    
    /** Updates the cumulative score function with a new onset detection function sample 
     * @param odfSample_p an onset detection function sample
     */
    void updateCumulativeScore (double odfSample_p);
	
    /** Predicts the next beat, based upon the internal program state */
    void predictBeat();
    
    /** Calculates the current tempo expressed as the beat period in detection function samples */
    void calculateTempo();
    
    /** Calculates an adaptive threshold which is used to remove low level energy from detection
     * function and emphasise peaks 
     * @param x_p a pointer to an array containing onset detection function samples
     * @param N_p the length of the array, x_p
     */
    void adaptiveThreshold (double* x_p, int N_p);
    
    /** Calculates the mean of values in an array between index locations [startIndex,endIndex]
     * @param array_p a pointer to an array that contains the values we wish to find the mean from
     * @param startIndex_p the start index from which we would like to calculate the mean
     * @param endIndex_p the final index to which we would like to calculate the mean
     * @returns the mean of the sub-section of the array
     */
    double calculateMeanOfArray (double* array_p, int startIndex_p, int endIndex_p);
    
    /** Normalises a given array
     * @param array_p a pointer to the array we wish to normalise
     * @param N_p the length of the array_p
     */
    void normaliseArray (double* array_p, int N_p);
    
    /** Calculates the balanced autocorrelation of the smoothed onset detection function
     * @param onsetDetectionFunction_p a pointer to an array containing the onset detection function
     */
    void calculateBalancedACF (double* onsetDetectionFunction_p);
    
    /** Calculates the output of the comb filter bank */
    void calculateOutputOfCombFilterBank();
	
    //=======================================================================

    /** An OnsetDetectionFunction instance for calculating onset detection functions */
    OnsetDetectionFunction m_Odf;
    
    //=======================================================================
	// buffers
    
    CircularBuffer m_OnsetDF;                 /**< to hold onset detection function */
    CircularBuffer m_CumulativeScore;         /**< to hold cumulative score */
    CircularBuffer m_SpectralPower;           // measure for volume, detect song begin/end and stop beat detection
    
    double m_ResampledOnsetDF[512];           /**< to hold resampled detection function */
    double m_Acf[512];                        /**<  to hold autocorrelation function */
    double m_WeightingVector[128];            /**<  to hold weighting vector */
    double m_CombFilterBankOutput[128];       /**<  to hold comb filter output */
    double m_TempoObservationVector[41];      /**<  to hold tempo version of comb filter output */
    double m_Delta[41];                       /**<  to hold final tempo candidate array */
    double m_PrevDelta[41];                   /**<  previous delta */
    double m_PrevDeltaFixed[41];              /**<  fixed tempo version of previous delta */
    double m_TempoTransitionMatrix[41][41];   /**<  tempo transition matrix */
    
	//=======================================================================
    // parameters
    
    double m_Tightness;                       /**< the tightness of the weighting used to calculate cumulative score */
    double m_Alpha;                           /**< the mix between the current detection function sample and the cumulative score's "momentum" */
    double m_BeatPeriod;                      /**< the beat period, in detection function samples */
    double m_Tempo;                           /**< the tempo in beats per minute */
    double m_EstimatedTempo;                  /**< the current tempo estimation being used by the algorithm */
    double m_LatestCumulativeScoreValue;      /**< holds the latest value of the cumulative score function */
    double m_TempoToLagFactor;                /**< factor for converting between lag and tempo */
    int m_m0;                                 /**< indicates when the next point to predict the next beat is */
    int m_BeatCounter;                        /**< keeps track of when the next beat is - will be zero when the beat is due, and is set elsewhere in the algorithm to be positive once a beat prediction is made */
    int m_HopSize;                            /**< the hop size being used by the algorithm */
    int m_OnsetDFBufferSize;                  /**< the onset detection function buffer size */
    bool m_TempoFixed;                        /**< indicates whether the tempo should be fixed or not */
    bool m_BeatDueInFrame;                    /**< indicates whether a beat is due in the current frame */
    int m_FFTLengthForACFCalculation;         /**< the FFT length for the auto-correlation function calculation */
    
    double m_PowerAtLastBeat = 0.;

#ifdef USE_FFTW
    //fftw_plan acfForwardFFT;                /**< forward fftw plan for calculating auto-correlation function */
    //fftw_plan acfBackwardFFT;               /**< inverse fftw plan for calculating auto-correlation function */
    fftw_complex* m_ComplexIn;                /**< to hold complex fft values for input */
    fftw_complex* m_ComplexOut;               /**< to hold complex fft values for output */
#endif
    
#ifdef USE_KISS_FFT
    kiss_fft_cfg m_CfgForwards;               /**< Kiss FFT configuration */
    kiss_fft_cfg m_CfgBackwards;              /**< Kiss FFT configuration */
    kiss_fft_cpx* m_FFTIn;                    /**< FFT input samples, in complex form */
    kiss_fft_cpx* m_FFTOut;                   /**< FFT output samples, in complex form */
#endif

};

#endif
