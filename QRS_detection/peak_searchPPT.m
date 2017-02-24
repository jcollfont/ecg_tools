function [peaks z] = peak_searchPPT(signal, Nw, tha, tav, initMax, searchT, idleT)
%% HELP:
%		[peaks z] = peak_searchPPT(signal, Nw, tha, tav, initMax, searchT, idleT)
%
%			This function searches for peaks in each of the dimensions the given input data.
%			The algorithm implemented is the presented in:
%				W Zong, T Heldt, GB Moody, RG Mark. "An Open-source
%				Algorithm to Detect Onset of Arterial Blood Pressure Pulses".
%				Computers in Cardiology, 2003.
%
%		INPUT:
%			- signal - <L,N>double - input signal.
%			- Nw - int - window size.
%			- tha - double - threshold a.
%			- tav - double - number of previous peaks used for max calculation.
%			- initMax - double - initial maximum to use before the
%			averaging.
%			- searchT - int - number of samples to search a maximum.
%			- idleT - int - idle time after detection;
%
%		OUTPUT:
%			- peaks - <L,NP>int - position of the peaks.
%
%		PROCESS:
%
%		DEPENDENCES:
%
%		AUTHOR:
%			
%
%
	
	%% DEFINE
		[L N ] = size(signal);
		peaks = zeros(L,0);
		
	%% for all the dimensions
	for l = 1:L
		
		%% prefilter signal
			
		
		%% calculate dy
			dy = signal(l,2:end) - signal(l,1:end-1);
			
		%% calculate du
			mask = zeros(size(dy));
			mask(dy>0) = 1;
			du = dy.*mask;

		%% Slope sum function
			
			% sort signal ina  matrix S
			S = zeros(Nw,N-1);
			for ii = 1:Nw
				S(ii,:) = [ zeros(1,ii-1) du(1:N - (ii)) ];
			end
			
			% calculate slope
			z = sum(S,1);
			
		%% Decision making
			maxim = initMax;
			start = 1;
			nMax = 1;

			while true

				cross = find( z(start:min(N-1,start + searchT)) >= tha*maxim );
				it = 1;
				tmpSearch = searchT;
				while numel(cross) == 0
					cross = find( z(start:min(N-1,start+ tmpSearch)) >= (1-0.05*it)*tha*maxim );
					if  (((1-0.05*it)*tha) < 0.1) || (it > 20) && (numel(cross) == 0)
						tmpSearch = tmpSearch + searchT;
						cross = [];
						it = 0;
					end
					it = it +1;
				end
				
				% search peak within the search time
				[val(nMax) m(nMax)] = max( z( cross(1) + start : min(N-1,cross(1)+ searchT + start) ) );
				
				
				peaks(l,nMax) = m(nMax) + start + cross(1);
% 				
% 				plot(z);hold on;plot(cross(1) + start:min(N-1,cross(1)+ searchT + start),z(cross(1) + start:min(N-1,cross(1)+ searchT + start)),'r'); plot(peaks(l,1:nMax),z(peaks(l,1:nMax)),'go');hold off;
% 				pause;
				
				% update maximum
				if nMax > tav
					maxim = mean(val(nMax-tav:nMax));
				end
				
				start = peaks(l,nMax) + idleT;
				if start > N - idleT
					break;
				end
				nMax = nMax + 1;
			end
			
		%% correct the peaks
		plot(z);
	end
end