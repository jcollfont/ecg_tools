function [filt_ECG, baseline, refPoints, refIntervals] = baselineCorrection_Splines_auto(ECG, margin, refSize, refract)
%% HELP:
%		[filt_ECG] = baselineCorrection_Splines(ECG)
%			This function applies baseline correction to the introduced ECG
%			signal with baseline correction.
%			The algorithm is semiautomatic. It requires the user to select
%			the 0 intervals on the rms of the input signals.
%
%		INPUT:
%			- ECG - <L,N>double - input ECG signal. It is structured in a
%			matrix of L leads by N time instances.
%			- margin - int - number of samples before the R wave where to
%			select the knot.
%			- refSize - int - size of the reference interval.
%			- refract - int - number of samples during which no beat should
%			be detected. (estimated length of beat)
%
%		OUTPUT:
%			- filt_ECG - <L,N>double - filtered ECG signal. Same structure
%			as the input.
%			- baseline - <1,L>cell - each cell contains:
%									<1,N>double - subtracted baseline.
%			- refPoints - <1,L>cell - each cell contains:
%									<2,NR>double - position of the baseline
%									reference points.
%			- refIntervals - <2,NR>cell - each cell contains the start or
%			the end of the selected interval to do calculate the baseline.
%
%		PROCESS:
%			- calculate rms of the signal.
%			- user selection of the reference ('0') intervals.
%			- for each lead:
%				- average the intervals to find the reference points.
%				- calculate the spline that passes through all ref points.
%				- substract estimated baseline.
%			- plot resulting rms.
%
%		DEPENDENCES:
%			- sepectPoints_onSignal.m
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%
%

	%% DEFINE
		[L N] = size(ECG);
		T = 5000;
		filt_ECG = ECG;
        baseline = cell(1,L);
        refPoints = cell(1,L);
		
		% check existance of t
		if ~exist('t')
			t = 1:1:N;
		end
		
		if ~exist('refract')
			refract = 100;
		end
		
	%% fix reference point (remove mean accross lead space)
		if L >1
				e = ones(1,L);
			ECG = ECG - repmat(1/L*e*ECG,L,1);
		end
	
	%% calcualte rms of the signal
		rmsECG = sqrt( sum(ECG.^2,1) );
		
	%% User selection of the intervals
		[Rwave] = rWave_Detect_Wavelet3(ECG,0,refract);
		
		nInt = numel(Rwave);
		refIntervals = cell(1,2);
		refIntervals{1} = zeros(1,nInt);
		refIntervals{2} = zeros(1,nInt);
		falseStart = false;
		for r = 1:nInt
			
			tmp = Rwave(r) - margin;
			if tmp >= 1
				refIntervals{1}(r) = tmp;
				refIntervals{2}(r) = tmp + refSize;
			else
				falseStart = true;
			end
			
		end
		
		if falseStart
			refIntervals{1} = refIntervals{1}(2:end);
			refIntervals{2} = refIntervals{2}(2:end); 
		end
			
	%% For each lead
		for l = 1:L
			%% Average intercals to find reference points
				NR = numel(refIntervals{1});
				refPoints{l} = zeros(2,NR);
				for nr = 1:NR
					refPoints{l}(1,nr) = mean( ECG(l, refIntervals{1}(nr):refIntervals{2}(nr)) );
					refPoints{l}(2,nr) = (refIntervals{2}(nr) + refIntervals{1}(nr) )/ 2;
				end
				
				%% fit splines to reference points to obtain baseline
					baseline{l} = spline(refPoints{l}(2,:),refPoints{l}(1,:), t);
				
				%% substract baseline from the original signal
					filt_ECG(l,:) = ECG(l,:) - baseline{l};
					
		end%end of all leads
		
		
	% plot resulting rms
% 		figure; plot(ECG(l,:)');plot(baseline{l} ,'r');
% 		rmsECG = sqrt( sum(filt_ECG.^2,1) );
% 		figure;title('filtered ECG rms');plot(rmsECG);
		
end% end of function