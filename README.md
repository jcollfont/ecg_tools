# ECG_TOOLS

This repo contains various tools for ECG signal processing used by Jaume Coll-Font.

The contents can be found in the following groups:
1. Beat segmentation and alignment
2. Denoising
3. QRS detection

The three blocks depend on each other, I do not recommend to separate them.

I have tried to describe each file in the headers, but this is mostly "research code", which means that, at times, might be incomplete or poorly described.

This code has been developed by Jaume Coll-Font.
If you find issues with it, please contact me at jcollfont ADD gmail.com


## Beat segmentation and alignment
This folder contains a series of functions to perform beat segmentation. Some of them are very specific to some datasets, but the general ones are:
* **segmentHeartBeats_bySections_fixedLength:** Segments all the heartbeats of a signal by selecting a window around the Q-peak (manually or automatically selected)
* **segmentTMPs:** Segments beats in Transmembrane Potentials signals. Uses the assumption that the depolarization wave is markedly steeper than anything else.
	* **segmentTMPs_SameCut:** applies the same approach but also selects the corresponding beats in extracellular potentials and body surface potentials (useful for datasets generated with ECGSIM).
* **segment_Twave_paramBased:** uses the parameters defined in an ECGSIM case to do the beat segmentation.
* **alignTwaves:** aligns a series of segmented T-waves based on a desired metric.


## Denoising
This folder contains a series of functions to reduce noise in an ECG signal. Some of them are very specific to some datasets, but the general ones are:
* **baselineCorrection_Splines_auto:** eliminates the baseline wander present in the ECG recordings.

## QRS detection
This folder contains a series of functions to detect the QRS. Some of them are very specific to some datasets, but the general ones are:
* **rWave_Detect_Wavelet3:** uses a wavelet based approach to detect the peak of QRS. The same approach could be used to detect other fiducials.
