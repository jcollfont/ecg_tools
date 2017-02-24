function [Tseg qrs_peaks cuts oSeg] = segmentHeartBeats_bySections_fixedLength(ECG, startT, lengthT, mode, doPlot, preMarkedPoints)
%% HELP:
%		[Tseg qrs_peaks cuts oSeg] = segmentHeartBeats_bySections_fixedLength(ECG, startT, lengthT, mode, doPlot)
%			This algorithm segments the heartbeats manually with a fixed
%			length.
%			To do so, the user must input the length of all Twaves and also
%			select the start of every T-wave.
%			The usage is left click on the vertical of the desired time
%			instance and right click to go to the next section/finish.
%
%		INPUT:
%			- ECG - <1,NS>cell - each cell contains a signal to segment. The first will be taken as 
%								as reference to define where to segment
%								from. Each signal has de form:
%								- <L,N>double - raw ECG signals. L are the number of
%								leads and N the total of samples.
%			- startT - int - average number of samples that takes the T
%			wave to start after the QRS.
%			- lengthT - int - average length of the Twaves.
%			- mode - int - desired segmentation mode:
%							0 - manual; the user manually selects the start
%							of the Twave.
%							1 - automatic; the algorithm selects the peak
%							of the QRS complex as reference. The selection
%							is made using wavelets.
%							2 - use the preselected qrs peaks from preMarkedPoints
%			- doPlot - int - if different than 0, plot the corresponding
%			lead.
%			- preMarkedPoints - <1,NB>int - (OPTIONAL) vector with the preselected qrs_peaks.
%
%		OUTPUT:
%			- Tseg -<1,NB>cell - each cell contains the Twave of each
%					heartbeat in format:
%					-<L,N>double - where L is the number of leads and N the
%					number of samples. (N must be equal for all beats)
%			- qrs_peaks - <1,NB>int - position of the detected QRS peaks.
%			- cuts - <2,NB>int - starting and end points of each Twave.
%
%		PROCESS:
%
%		DEPENDENCES:
%			- selectPointds_onSignal.m
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%
%
%
	%% DEFINE
		if iscell(ECG)	% this solution is to keep compatibility with previous versions
			NS = numel(ECG);
			oSignals = ECG(2:end);
			ECG = ECG{1};
			oSeg = cell(1,NS-1);
			multSig = true;
		else
			NS = 0;
			oSeg = {};
			oSignals = {};
			multSig = false;
		end
	
		figT = 5000;		
		[L, N] = size(ECG);
		
		if ~exist('doPlot')
			doPlot = false;
		end
	
	switch(mode)
		case 0
			%% calculate rms of the signal
				fprintf('\t Calculating rms over lead space of the ECG for manual segmentation.\n');
				rmsECG = sqrt(sum(ECG.^2,1));

			%% select start of the Twaves
				fprintf('\t Select the intervals where the signal should be 0V (type 1 is start, type 2 is end of interval).\n');
				[segPoints] = selectPoints_onSignal(rmsECG, figT, 1);
				
		case 1
% 			[QRS_points] = qrsDetect_Wavelet(ECG, [1:N]);
			%% calculate rms of the signal
				fprintf('\t Calculating rms over lead space of the ECG for manual segmentation.\n');
				rmsECG = sqrt(sum(ECG.^2,1));
			
			%% search QRS peaks
				[rwave] =   rWave_Detect_Wavelet3(ECG,false,lengthT+startT);
				segPoints{1} = rwave;
				
		case 2
			segPoints{1} = preMarkedPoints;
			
	end
	

	%% apply segmentation
		fprintf('\t Applying segmentation.\n');
		numSeg = numel(segPoints{1});
		if segPoints{1}(end) +startT + lengthT > N
			numSeg = numSeg -1;
		end
		
		if multSig	% case we want to segment multiple signals...
			for ii = 2:NS
				oSeg{ii-1} = cell(1,numSeg);
			end
			
			Tseg = cell(1,numSeg);
			cuts = zeros(2,numSeg);
			for ii = 1:numSeg
				Tseg{ii} = ECG(:,segPoints{1}(ii)+startT:segPoints{1}(ii)+startT+lengthT-1);
				cuts(:,ii) = [segPoints{1}(ii)+startT ; segPoints{1}(ii)+startT+lengthT-1];
				for jj = 2:NS
					oSeg{jj-1}{ii} = oSignals{jj-1}(:,segPoints{1}(ii)+startT:segPoints{1}(ii)+startT+lengthT-1);
				end
			end
			
		else % just one signal
			Tseg = cell(1,numSeg);
			cuts = zeros(2,numSeg);
			for ii = 1:numSeg
				Tseg{ii} = ECG(:,segPoints{1}(ii)+startT:segPoints{1}(ii)+startT+lengthT-1);
				cuts(:,ii) = [segPoints{1}(ii)+startT ; segPoints{1}(ii)+startT+lengthT-1];
			end
		end
		
		
		
		qrs_peaks = segPoints{1};
		
		
	%% plot
		if doPlot
			l= doPlot;
			figure;hold on;
			plot(ECG(l,:),'k--')
			for ii = 1:numSeg
				t = cuts(1,ii):cuts(2,ii);
				plot(t,Tseg{ii}(l,:),'r')
			end
			plot(qrs_peaks,ECG(l,qrs_peaks),'go');
			hold off;
		end

end