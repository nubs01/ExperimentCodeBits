--------------------------------------------------------------------------------
                   ___                __           _     
                  /   |  ____  ____ _/ /_  _______(_)____
                 / /| | / __ \/ __ `/ / / / / ___/ / ___/
                / ___ |/ / / / /_/ / / /_/ (__  ) (__  ) 
               /_/  |_/_/ /_/\__,_/_/\__, /____/_/____/  
                                    /____/               
              _       __           __   ______             
             | |     / /___  _____/ /__/ __/ /___ _      __
             | | /| / / __ \/ ___/ //_/ /_/ / __ \ | /| / /
             | |/ |/ / /_/ / /  / ,< / __/ / /_/ / |/ |/ / 
             |__/|__/\____/_/  /_/|_/_/ /_/\____/|__/|__/  
                                                           
--------------------------------------------------------------------------------

# Indiviual Animal Analysis

## Add to preprocess
- include gender in rntask file
- function to update animal direction if LEDs are on sides rather than front-back
- make state struct for each day containg state params and start & end times of behavioral states
- detect sleep spindles and make struct like ripple struct
- check data validity
    - Identify artifacts (make data structure with artifact times)
    - Check position tracking time vs LFP time (before preprocessing)
    - Check diode numbers (before preprocessing) - automatically generate diodenums for preprocessing
- make plots:  
    - position tracking overview (colored for state)
    - overall spectrum
    - overall spectrogram with colored lines above indicating state for every day/riptet w/ epoch dividers
    - average state spectra (bootstrapped SEM) - raw & normalized for every epoch/riptet
    - average state spectra (bootstrapped SEM) - raw & normalized over all riptets for each epoch
    - subplots with EMG + envelope, cleaned velocity, behavioral state (like a hypnogram)
    - subplots for each power band: for each state transition (moving -> sleep, etc) order episodes by length and plot on time (-5s to Ns) vs episode length with lines colored for band power



# Gather data to analyze

## Animal Metadata
-  Animal IDs
-  Genotypes
-  Genders
-  Data directories

## Analysis Metadata
- Analysis Directory
- Window Size
- Spectral Window Size
- Behavioral State Parameters 
- Power Bands to analyze (band-freqs & band-names) - delta, theta, sigma, gamma, 
- Normalization Function (@(x,y,z) (x-y)/z) :  (band_power_in_window - median_power_in_epoch)/stdev_of_power_in_epoch
- Normalization Band & Method (total vs mean vs median) : 1-100 Hz

## Collect Unified Dataset with metrics for each window for each animal/day/epoch/tet
Return Windowed Data:
    - Low Freq Power Spectrum from EEG
    - High Freq Power Spectrum from EEG-Ref
    - # of ripples in segment
    - Behavioral State
    - mean (or median), st.dev & SEM(bootstrapped) for normalization band computed over whole epoch (EEG)
    - mean (or median), st.dev & SEM(bootstrapped) for normalization band computed over whole epoch (EEG-REF)
    - Raw band powers (delta,theta,gamma) for segment
    - Max velocity in segment
    - Animal, Day, Epoch, Tetrode, segment-time

Return Epoch Data:
    - Spectrogram Low Freq from EEG
    - Spectrogram High Freq from EEG-REF
    - Behavioral State Matrix
    - Ripple times
    - Artifact times
    - Total Distance
    - Phase-locking of cells (delta, theta, sigma, and gamma) - matrix of spike time and phase of every band
    - Spindle times
    - Animal, Day Epoch, Tetrode

Return data for each segment:
    - Animal, Day, Epoch, Tetrode, Segment_Time
    - lo_band_powers - EEG
    - hi_band_powers - EEG-REF
    - normalization power & SEM - EEG
    - normalization power & SEM - EEG-REF
    - ripples times in segment
    - artifact times in segment
    - behavioral state

# Plotting and Comparisons
Color Df1 with blue scale and WT with red scale
## Mean Power Spectrum for each animal - state independent -- bootstrapped SEM
## Mean Power Spectrum for each animal per state (sleeping, resting, moving, not-moving)
## Mean Power Spectrum for each group - state independent
## Mean Power Spectrum for each group per state (sleeping, resting, moving, not-moving)
## Bootstrapped means for all animals per state & power band
## Boostrapped means for each group per state & power band
## Fold change plot for state transitions (resting -> sleeping & resting -> moving) for each power band
## Bar plot of baseline power for each animal for each band (one plot)
## Bar plot of baseline power for each group for each band (one plot)
## Histogram of MUA phase-locking (# of spikes per phase bin) - for each power band


## Sal-Ket Analysis
- For each tetrode/day grab delta band power traces for sleep, sal & ket
- Normalize all to mean and std of sleep epoch 
- Make plot and then store normalized traces to average over all tetrodes

# Using spike data -- for each animal and for each group

## Phase-locking to delta, theta and gamma bands (single unit)
## Spike-phase histogram for delta, theta and gamma bands (multi-unit)
## Sike-phase histograms (single unit)
## Theta phase precession 
## 
