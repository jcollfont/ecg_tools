function  [ins_freq] = instant_freq_basedOnQRS(ECG, t, window)
%% HELP:
%		[ins_freq] = instant_freq_basedOnQRS(ECG, t, window)
%			This algorithm detects the instant frequency of the ECG signal.
%			It uses the qrs_detector over an user-defined window to
%			estimate the instant frequency.
%
%		INPUT:
%			- ECG - <L,N>double - input ECG.
%			- t - <1,N>double - sampling times of the ECG.
%			- window - int - size of the desired window.
%
%		OUTPUTS:
%			- ins_freq - <1,N>double - instant frequency.
%
%		PROCESS:
%
%		DEPENDENCES:
%			- qrsDetect_wavelet.m
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%
%
	
	%% DEFINE
		[L, N] = size(ECG);
		
		
	%% detect QRS
		[QRS_peak freq] = qrsDetect_Wavelet(ECG, t);
		
		
	%% for all time instances
	for ii = 1:N
		
		%% window
			if floor(window/2) > ii -1
				
				startT = 1;
				endT = ii + floor(window/2);
			elseif floor(window/2) > N - ii
				
				startT = ii - floor(window/2) +1;
				endT = N;
			else
				
				startT = ii - floor(window/2) +1;
				endT = ii + floor(window/2);
			end
				
		%% find corresponding peaks
			IX = find( (QRS_peak{1} > startT) & (QRS_peak{1} <endT) );
			
		%% measure frequency
			if numel(IX) < 2
				freq = 0;
			else
				freq =  mean( 1./(t(QRS_peak{1}(IX(2:end))) - t(QRS_peak{1}(IX(1:end-1))) ) );
			end
			
			ins_freq(ii) = mean(freq(:,1),1);
			
	end
	
end