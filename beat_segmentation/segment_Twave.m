function [] = segment_Twave(tmp, fhPot, fbPot, window)
%% HELP:
%		segment_Twave:
%			This function automatically segments the Twave based on the
%			body surface potentials and extends the segmentation points to
%			heart potentials and TMPs.
%			This function uses the rms of the torso potentials to fins the
%			Twave.
%
%		INPUT:
%			- tmp - <M,N>double - matrix containing the tmp series of N
%			samples for the M points on the heart.
%			- fhPot - <M,N>double - matrix containing the heart potentials series of N
%			samples for the M points on the heart.
%			- fbPot - <L,N>double - matrix containing the torso potentials series of N
%			samples for the L points on the torso.
%
%		OUTPUT:
%			- stmp - <1,M>cell - each cell contains a beat of the tmp.
%			- sfhPot - <1,M>cell - each cell contains a beat of the heart potentials.
%			- sfbPot - <1,M>cell - each cell contains a beat of the torso potentials.
%			- cuts - <1,M+2>int - separation indices.
%
%		PROCESS:
%
%		DEPENDENCES:
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%
%
	
	%% DEFINE
		rms = sum(fbPot.^2,1);
		pThreshold = 0.4;
	
	%% calculate derivative
		derRms = rms(2:end) - rms(1:end-1);
		
		derRmsMtrx = repmat(derRms,window,1);
		mask = ones(size(derRmsMtrx));
		for w = 2:window
			mask(w,1:w-1) = 0;
			derRmsMtrx(w,:) = [zeros(1,w-1) derRms(1:end-w+1) ];
		end
			
	%% sum over positive
		posSumRms = sum((derRmsMtrx > 0).*derRmsMtrx,1);
		
	%% sum over negative
		negSumRms = abs(sum((derRmsMtrx < 0).*derRmsMtrx,1));
		
	%% find the peaks
		maxNeg = max(negSumRms);
		threshold = maxNeg*.4;
		IX = find(negSumRms == threshold);
		
	%% Find the onset
		
		
	%% find min/max
		figure;
		hold on;
		plot(rms,'b');
% 		plot(derRms,'m');
% 		plot(posSumRms,'r');
		plot(negSumRms,'g');
		hold off;
	
end