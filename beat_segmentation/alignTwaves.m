function [delays, t, alSignals, delaysLEADS] = alignTwaves(Tseg, type, signals, delVal)
%%HELP:
%	[delays t alSignals] = alignTwaves(Tseg, type, signals)
%		This function is thought to get a set of Twaves segmented from ECG
%		recordings and align them together.
%		The method used for this alignment is maximizing the joint cross
%		correlation of all leads in the recordings for each heartbeat.
%
%	INPUT:
%		- Tseg - cell<1,NB> - each cell represents a different heartbeat
%		and contains:
%				-<L,Nt>double - Twave signal. Where L is the number of leads and Nt the
%					number of samples. (Nt must be equal for all beats)
%		- type - string - [optional] type of alignment to apply. The default is
%		averaging. There are 3 types:
%				- mean: calculates the delay based on max xcorr of the mean
%				across leads.
%				- l1: calculates the delay based on max xcorr of the l1 norm
%				across leads.
%				- squared: calculates the delay based on max xcorr of the
%				l2 norm across leads.
%				- lead #: calculates the delay based on max xcorr of a single
%				lead.
%				- none: 0 delay applied across all beats.
%				- alllead: each lead is aligned separately according to its
%				max xcorr.
%				- max_l1: alings the maximum of the Twave l1 norm across leads.
%				- max_squared: alings the maximum of the Twave l2 norm across leads.
%				- fixed: it fixes the delay to a desired value in delVal.
%				- max_allead: picks the alignment based on the maximum of
%				each lead individually.
%		- signals - <1,NS>cell - each cell contains a set of signals the
%		user wants to align given the delays from Tseg. Each signal must
%		have the same structure as Tseg.
%		- delVal - int - (OPTIONAL) it is only necessary when type='fixed'.
%		It determines the fixed delay to apply to the T-waves.
%
%	OUTPUT:
%		- delays - <NB,L>double - delays of all beats with respect to the
%		first heartbeat.
%		- t - <1,NB>cell - each cell represetns a heartbeat and contains:
%				- <L,Nt>int - time stamps relative to the first heartbeat. (all in sample number)
%		- alSignals - <1,NS>cell - each cell contains a set of aligned signals the
%		user provided.
%		- delaysLEADS - <NB,L>double - if 'allleads' or 'max_allleads' selected, returns the
%		individual delay at each body surface lead.
%
%	PROCESS:
%		- cross correlate the first signal with all the others.
%		- find the maximum and the corresponding delay.
%		- create time vectors.
%
%	DEPENDENCES:
%
%	AUTHOR:
%		Jaume Coll-Font <jcollfont@gmail.com>
%
%

	%% Define
		[NB] = numel(Tseg);
		[L Nt] = size(Tseg{1});
		delays = zeros(NB,L);
		t = cell(1,NB);
		avrg = cell(1,NB);
		delaysLEADS = zeros(NB,L);
		
		if nargin == 1
			typeIX = 1;
				
		elseif nargin >= 2
				
			if strcmp(type,'mean')
				typeIX = 1;
			elseif strcmp(type,'squared')||strcmp(type,'l2')
				typeIX = 2;
			elseif strcmp(type,'none')
				typeIX = 4;
			elseif strcmp(type,'alllead')
				typeIX = 5;
			elseif strcmp(type,'max_l1')
				typeIX = 6;
			elseif strcmp(type,'max_l2')
				typeIX = 7;
			elseif strcmp(type,'l1')
				typeIX = 8;
			elseif strcmp(type,'fixed')
				typeIX = 9;
			elseif strcmp(type,'max_allead')
				typeIX = 10;
			else
				typeIX = 3;
				lead = str2num(type);
			end
		end
			
			
		
		
	%% For all heartbeats
		if (typeIX == 1)||(typeIX ==2)||(typeIX == 3)||(typeIX==8)
			for ii = 1:NB

				%% take average along the spatial dimension
					switch typeIX
						case 1
							avrgT{ii} = mean(Tseg{ii},1);
						case 2
							avrgT{ii} = sqrt(sum(Tseg{ii}.^2,1));
						case 3
							avrgT{ii} = Tseg{ii}(lead,:);
						case 8
							avrgT{ii} = mean(abs(Tseg{ii}),1);
					end

				%% align beats wrt the first beat
					[corre lag] = xcorr(avrgT{1},avrgT{ii});

				%% find the maximum x-correlation
					[sink ix] = max(abs(corre));
					delays(ii,:) = repmat(lag(ix),1,L);

			end
			delaysLEADS  = delays;
			
		elseif (typeIX == 5)
			for jj = 1:NB
				for l = 1:L
					[c lag] = xcorr( Tseg{1}(l,:),Tseg{jj}(l,:) );
					[sink IX] = max(abs(c));
					delaysLEADS(jj,l) = lag(IX);
				end
			end
			
			delays = repmat( round(mean(delaysLEADS,2)) ,1,L);
			
		elseif ( typeIX == 10 )
			[sink refIX] = max(abs(Tseg{1}),[],2);
			delaysLEADS(1,:) = 0;
			for jj = 1:NB
				[sink IX] = max(abs(Tseg{jj}),[],2);
				delaysLEADS(jj,:) = refIX - IX;
			end

			delays = repmat( round(mean(delaysLEADS,2)) ,1,L);
			
		elseif (typeIX == 6)||(typeIX == 7)
							
			% for all heartbeats
			for ii = 1:NB
				
				% take average along the spatial dimension
					switch typeIX
						case 6
							avrgT{ii} = mean(abs(Tseg{ii}),1);
						case 7
							avrgT{ii} = sqrt(sum(Tseg{ii}.^2,1));
					end
					
				% search for Twaves peaks
					[sink maxIX] = max(abs(avrgT{ii}));
					[sink maxRef] = max(abs(avrgT{1}));
					
				% set delay
					delays(ii,:) = repmat( maxRef - maxIX , 1, L);
			end
		elseif (typeIX == 9)
			
			delays(2:NB,:) = delVal;
			delaysLEADS  = delays;
			
		end
		
	if (typeIX ~= 5)&&(typeIX ~= 10)
		%% create time vectors
			minDel = min(min(delays));
			N = Nt + max(max(delays)) - minDel;
			for ii = 1:NB
				for ll =1:L
					% adjust individual delay wrt heartbeat 1
					t{ii} = repmat(1:N,L,1) + minDel - repmat(delays(ii,:)',1,N);
				end
			end

		%% align signals
			alSignals = {};
			if nargin >= 3

				NS = numel(signals);
				alSignals = cell(1,NS);

				for ii = 1:NS
					[L sink] = size(signals{ii}{1});
					alSignals{ii} = zeros(L,N,NB);
					for b = 1:NB
						for l = 1:L
							alSignals{ii}(l, (t{b}(1,:)>0)&(t{b}(1,:)<=Nt) ,b ) = signals{ii}{b}(l,:);
						end
					end
				end
			end

	else
		%% create time vectors
			minDel = min(min(delaysLEADS));
			maxDel = max(max(delaysLEADS));
			N = Nt + maxDel - minDel;
			for ii = 1:NB
				for ll =1:L
					% adjust individual delay wrt heartbeat 1
					t{ii}(ll,:) = [0:N] + minDel + 1 - delaysLEADS(ii,ll);
				end
			end

		%% align signals
			alSignals = {};
			if nargin >= 3

				NS = numel(signals);
				alSignals = cell(1,NS);

				for ii = 1:NS
					[Ls sink] = size(signals{ii}{1});
					alSignals{ii} = zeros(Ls,N,NB);
					for b = 1:NB
						if Ls ~= L
							for l = 1:Ls
								alSignals{ii}(l, (mean(t{b},1)>0)&(mean(t{b},1)<=Nt) ,b ) = signals{ii}{b}(l,:);
							end
						else
							for l = 1:Ls
								alSignals{ii}(l, (t{b}(l,:)>0)&(t{b}(l,:)<=Nt) ,b ) = signals{ii}{b}(l,:);
							end
						end
					end
				end
			end
			
	end

end