function [stmp t] = segmentTMPs(tmp)
%% HELP:
%		segmentTMPs:
%			segments the given TMPs. To do that, the algorithm seaches for
%			the midd point between two consecutive heartbeats.
%			This algorithm asumes a clear resting potential for the TMPs
%			with similar amplitude values for all beats.
%			Also asumes that the signal starts at resting potential;
%			there's no cutted beat at the begining. And all the points on
%			the heart will have the same ammount of heartbeats
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
	
	%% For all Leads
	for l = 1:L
		
		%% Search min and max value in the TMP
			maxV = max(tmp(l,:));
			minV = min(tmp(l,:));
			amplV = maxV - minV;
			
		%% Search for the beat searations
			tr = find( ((1-decayV)*(amplV+minV)+minV) - tmp(l,:) > 0 );
			jumpS = tr((tr(2:end) - tr(1:end-1)) > 1);
			jumpE = tr(find( (tr(2:end) - tr(1:end-1)) >1 )+1);
			cuts = ceil(( jumpS(2:end) + jumpE(1:end-1) )/2);
			
		
		%% Sort heart beats
			M = numel(cuts);
			tempSTMP{l} = cell(1,M +1);
			tempT{l} = cell(1,M+1);
			
			tempSTMP{l}{1} = tmp(l,1:cuts(1));
			tempT{l}{1} = 1:cuts(1);
			for m = 2:M-1
				tempSTMP{l}{m} = tmp(l,cuts(m)+1:cuts(m+1));
				tempT{l}{m} = (cuts(m)+1):cuts(m+1);
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
	
end