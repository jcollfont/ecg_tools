function [Tdom] = dominant_Twave(Tseg)
%% HELP:
%		dominant_Twave.m
%			This function evaluates the dominant Twave as descived by 
%			van Ooseterom in:
%			van Oosterom. "The Dominant T wave and its significance",
%			Journal of Cardiovascular Electrophysiology, Vol 14, No. 10
%			Supplement, October 2003.
%
%		INPUT:
%			- Tseg -<1,NB>cell - each cell contains the Twave of each
%			heartbeat in format:
%					-<L,NT>double - where L is the number of leads and NT the
%					number of samples. (NT must be equal for all beats)
%
%		OUTPUT:
%			- Tdom -<1,NT>double - dominant Twave of NT samples.
%
%		PROCESS:
%
%		DEPENDENCES:
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail..com>
%
%
	
	%% DEFINE
		NB = numel(Tseg);
		L = size(Tseg{1},1);
		
		% find maximum number of samples
		N = zeros(1,NB);
		for ii = 1:NB
		 [N(ii)] = size(Tseg{ii},2);
		end
		
		mN = max(N);
		
	%% create a matri with all the T waves staked together with the same number of samples (zero padding)
		Tmtrx = zeros(NB*L,mN);
		for ii =1:NB
			Tmtrx(L*(ii-1)+1:ii*L,:) = [Tseg{ii} zeros(L,mN - N(ii))];
		end
		
	%% svd the matrix
		[U S V] = svd(Tmtrx);
		
		Tdom = V(:,1)';
		
	
end