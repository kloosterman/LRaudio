This Readme file explains the structure of the single subject data recorded during the EEG experiment (spatial selectvie auditory attention paradigm).
The data are provided as single Matlab data files.

Reference	Woestmann, Alavash, Obleser, 2019 J Neurosci
Contact 	Mohsen Alavash, U Luebeck mohsen.alavash@uni-luebeck.de
Last update	10.10.2022
--------------------------
Each participant has completed six blocks of the task.
For each participant there is one sub_*.mat file. 
This file contains clean artefact-rejected EEG sensor data with single-trial behavioral measures and stimulus conditions attached to it.
The structure of the data follows the format used in Fieldtrip toolbox.
Each task block contains 96 trials (or less if trials were rejected due to poor data quality duriong preprocessing).
This gives around 6x96=576 (or less) trials per participant.

Details on preprocessing and data structure:
       	- Reference channel: 		online (TP9)
	- Band-pass filter: 		1-100 Hz
	- Smapling rate: 		250 Hz
       	- Epoch: 		-2s to 6s relative to spatial cue onset
       	- Side of auditory spatial cue: 		eegdata.trialinfo, 1st column
       	- Speaker arrangement: 		eegdata.trialinfo, 2nd column
       	- Accuracy and RT (relative to the stimulus offset): 		eegdata.trialinfo, 3rd and 4th columns
       	- Button press: 		eegdata.trialinfo, 5th column (65,70: high confidence) 		
  	- Tracking parameter (delta_f in semitones): 		eegdata.trialinfo, 6th column
	- Jitter (anticipation period in ms): 		eegdata.trialinfo, 7th column
	- Congruency of tone streams: 		eegdata.trialinfo, 8th column

For a few subjects (and blocks) it was not possible to match the trial information and EEG data due to technical issues. 
For these cases data are not available (NaNs).