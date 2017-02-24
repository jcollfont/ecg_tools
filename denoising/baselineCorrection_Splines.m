function [filt_ECG, baseline, refPoints, refIntervals] = baselineCorrection_Splines(ECG, t, refIntervals)
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
%			- t - <1,N>double - time stamp of the ECH samples. If not
%					available, it is assumed to be t = 1:1:N
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
		
	%% fix reference point (remove mean accross lead space)
		if L >1
				e = ones(1,L);
			ECG = ECG - repmat(1/L*e*ECG,L,1);
		end
	
	%% calcualte rms of the signal
		rmsECG = sqrt( sum(ECG.^2,1) );
		
	%% User selection of the intervals
		if ~exist('refIntervals')
			[refIntervals] = selectPoints_onSignal(rmsECG, T, 2);
		end
		fprintf('Select the intervals where the signal should be 0V (type 1 is start, type 2 is end of interval).\n');
		% check the intervals are closed
			if numel(refIntervals{1}) ~= numel(refIntervals{2})
				fprintf('Invalid selected intervals.\n');
				return;
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