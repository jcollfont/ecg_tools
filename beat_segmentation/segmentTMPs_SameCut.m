function [stmp t cuts sHPot sBPot] = segmentTMPs_SameCut(tmp, fhPot, fbPot)
%% HELP:
%		segmentTMPs:
%			segments the given TMPs. To do that, the algorithm seaches for
%			the midd point between two consecutive heartbeats.
%			In this version the cut point is equal for all the nodes on the
%			heart.
%			This algorithm asumes a clear resting potential for the TMPs
%			with similar amplitude values for all beats.
%			Also asumes that the signal starts at resting potential;
%			there's no cutted beat at the begining. And all the points on
%			the heart will have the same ammount of heartbeats.
%			Another assumption is that the first depolariation of a beat
%			does not overlap with the last repolarization of the prvious
%			beat.
%
%		INPUT:
%			- tmp - <M,N>double - matrix containing the tmp series of N
%			samples for the M points on the heart.
%
%		OUTPUT:
%			- stmps - <L,M>cell - each cell contains the tmp of the
%			heartbeat corresponding to the l point on the heart at the mth
%			beat.
%			- t - <L,M>cell - each cell contains the time of the
%			heartbeat corresponding to the l point on the heart at the mth
%			beat.
%			- cuts - <1,M>int - cut indices used.
%
%		PROCESS:
%			For all Leads
%				Search min and max value in the TMP
%				Search for the beat separations
%				Sort heart beats.
%				
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%

	%% DEFINE
	[L N] = size(tmp);
	decayV = 0.9; %percentage of decay in amplitude
	tempSTMP = cell(1,L);
	tempT = cell(1,L);
	
	cutForwardPot = false;
	if nargin > 1
		cutForwardPot = true;
	end
	
	%% For all Leads
	for l = 1:L
		
		%% Search min and max value in the TMP
			maxV = max(tmp(l,:));
			minV = min(tmp(l,:));
			amplV = maxV - minV;
			
		%% Search for the beat searations
			tr = find( ((1-decayV)*(amplV+minV)+minV) - tmp(l,:) > 0 );
			jumpS(l,:) = tr((tr(2:end) - tr(1:end-1)) > 1);
			jumpE(l,:) = tr(find( (tr(2:end) - tr(1:end-1)) >1 )+1);
			
			
	end
	
	%% Find a common cutting point
		genJumpS = min(jumpS);
		genJumpE = max(jumpE);
		
		cuts = ceil(( genJumpS(2:end) + genJumpE(1:end-1) )/2);
		
		[s M] = size(cuts);
	
	%% For all Leads; CUT
	for l = 1:L
		
			%% Sort heart beats
			
			tempSTMP{l} = cell(1,M +1);
			tempT{l} = cell(1,M+1);
			
			tempSTMP{l}{1} = tmp(l,1:cuts(1));
			tempT{l}{1} = 1:cuts(1);
			for m = 1:M-1
				tempSTMP{l}{m+1} = tmp(l,cuts(m)+1:cuts(m+1));
				tempT{l}{m+1} = (cuts(m)+1):cuts(m+1);
			end
			tempSTMP{l}{M+1} = tmp(l,cuts(M)+1:end);
			tempT{l}{M+1} = cuts(M)+1:numel(tmp(l,:));
	end
	
	stmp = cell(L,M+1);
	t = cell(L,M+1);
	%% join all leads
	for l = 1:L
		for m = 1:M+1
			stmp{l,m} = tempSTMP{l}{m};
			t{l,m} = tempT{l}{m};
		end
	end	
	
	
	
	%% Do the heart and body surface segmentation if asked
	sHPot = {};
	sBPot = {};
	if cutForwardPot
		
		% Heart surface
		[L N] = size(fhPot);
		sHPot = cell(L,M+1);
		for l = 1:L
			sHPot{l,1} = fhPot(l,1:cuts(1));
			for m = 1:M-1	
				sHPot{l,m+1} = fhPot(l,cuts(m)+1:cuts(m+1));
			end
			sHPot{l,M+1} = fhPot(l,cuts(M)+1:end);
		end
		
		% Body surface
		[L N] = size(fbPot);
		sBPot = cell(L,M+1);
		for l = 1:L
			sBPot{l,1} = fbPot(l,1:cuts(1));
			for m = 1:M-1	
				sBPot{l,m+1} = fbPot(l,cuts(m)+1:cuts(m+1));
			end
			sBPot{l,M+1} = fbPot(l,cuts(M)+1:end);
		end
	end
	
	

	
end