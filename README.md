# Filtering GUI

This is a MATLAB based GUI aimed at helping people with little to no programming knowledge perform filtering and ridge extraction on wavelets.

This program can read a signal, calculate it's wavelet transform and plot it's average amplitude/power. Then once can select the frequency band to be used for the extraction/filtering and also look at the fourier spectrum of the signal. 

This release currently needs to be run on MATLAB. An installer which does not need MATLAB will be put up soon.

## Getting Started
Download all the files in this repository along with all the ones in the auxilaryFunctions repository. Add both to path in MATLAB and simply run Filtering.m

## Instructions
#### Importing a signal
Click on File to load your signal in .csv or .mat format
A good signal for testing would be the bloodflow signal which can be found in the signals folder. Orientation for all signals found in this repository is rowwise

#### Selecting Data
- Use the zoom tool to select your data. Once you have chosen it click on Refresh(CLicking on refresh or deselcting the tool tip is necessary to update values)
- If you would like to type in the data limits. Simply type it into the box labelled Xlim

#### Refreshing a plot
Simply click on the signal fromt the signal list to refresh the wavelet plot.

#### Wavelet Transform Options
- Sampling Frequency: Rate at which the signal was sampled
- Max Frequency: Maximal frequency for which to calculate Wavelet Transform. Default value is sampling_frequency/2
- Min Frequency: Minimal frequency for which to calculate Wavelet Transform. Default value is the minimal frequency for which at least one WT coefficient is determined
- Central Frequency: Central frequency of the wavelet
- Under Sampling Rate: The number indicates every nth point that is used for plotting. WARNING: Going too high can lead to erroneous results
- Amplitude/Power: Plots the amplitude or the power plot of the wavelet transform

#### Band Marking
- Mark region: A cross hair will show up, now clicking anywhere on the wavelet transform plot will mark a line. Clicking again marks the second limit. There is no particular order in which they need to be marked. After satisfactory marking of the region, click on add region to add it to the list.
- Interval List: To delete a region from the list saimply highlight the region and press the delete key on your keyboard.
- Frequency boxes: One can also simply type in the limits to be precise and click enter after typing into the box. This will show the region on the wavelet transform. Use the add region button to add the frequency pair to the list.

#### Advanced Options
If you don't know what you're doing. DO NOT TOUCH THIS
- Wavelet Type: The type of wavelet used for the wavelet transform
- Preprocess: Subtract the 3rd order polynomial fit and then bandpass the signal in the band of interest [fmin,fmax]
- CutEdges: Determine if WT coefficients should be set to NaNs out of the cone of influence
- Extraction type: Choose between the butterworth filter or ridge extraction from the NMD tool box

#### Buttons
- Calculate Transform: Calculates and plots the wavelet transform
- Filter Signal:Filters the signal and 
- 2D Plot: Rotates the view point of the surf plot to the xy-plane 
- Intervals: Do not use. Under construction.

#### Plotting
When viewing the bands and the fourier plot. One can select multiple bands which will superimpose on the same plot. 
To do this hold down control and click on the desired bands. Shift click also works.


