function [Rwave] = rWave_Detect_Wavelet3(x,optPlot,REFRAC)
%%	HELP:
%		[QRS_peak freq] = rWave_Detect_Wavelet3(x,t)
%			This function detects the peak of the R wave in the input ECG.
%			It uses the rms of the input signal to calculate it with a
%			wavelet transform method based on:
%			Juan Pablo Martínez et al., "A Wavelet-Based ECG Delineator:
%			Evaluation on Standard Databases", IEEE transactions on
%			biomedical engineering, April 2004.
%
%		INPUT:
%			- x - <L,N>double - input ECG.
%			- optPlot - boolean - true if the user wants to plot the
%			results.
%
%		OUTPUT:
%			- Rwave - <Nqrs,1>int - position of the R-waves on the input
%			ECG.
%
%		PROCESS:
%			- Define
%			- define wavelet coeffitiens
%			- calculate rms of the signal
%			- calculate wavelet transform
%			- define thresholds
%			- search for R waves
%			- plot if asked for
%
%		DEPENDENCES:
%
%		AUTHOR:
%			Jaume Coll-Font <jcollfont@gmail.com>
%
%

	%% define
		SEP = 20; %max separation between jumps
		JJ = 10; % num of sampled scales
		[L, N] = size(x);
		QRS = cell(1,L);
		
		numBlocks = 5;
		
		if ~exist('optPlot')
			optPlot = false;
		end
		
		switch nargin
			case 1
				optPlot = false;
			case 2
				REFRAC = 100; % refractory period
		end
			
	
	%% Define wavelet coeffitiens
		h = [1 3 3 1 0]/8;
		g = [0 1 -1 0 0]*2;
	
	%% calculate rms of the signal
		x = sum(x.^2,1);
	
	%% calculate wavelet transform
		Wx = {x};Sx={x};for jj = 1:JJ; Wx{jj+1} = filter(g,1,Sx{jj}(1:2:end));Sx{jj+1} = filter(h,1,Sx{jj}(1:2:end));end;
			
		% correct subsampling in WT
		for ii = 1:JJ
			tmp = Wx{ii+1};
			for jj = 1:ii
				NZ = numel(tmp);
				tmp = reshape(repmat(tmp,2,1),1,2*NZ);
				tmp = [tmp(4:end) zeros(1,4-1)];
			end
			Wx{ii+1} = tmp;
		end
	
	%% define thresholds
% 		th = zeros(1,4);
% 		for ii =1:3
% 			th(ii) = sqrt( sum(Wx{ii+1}.^2)/N );
% 		end
% 			th(4) = .5*sqrt( sum(Wx{ii+1}.^2)/N );
			
		th = zeros(4,numBlocks);
		for jj = 1:numBlocks
			times = [1:ceil(N/numBlocks)] + (jj-1)*ceil(N/numBlocks);
			if jj == numBlocks
				times = ((jj-1)*ceil(N/numBlocks)+1):N;
			end
			for ii =1:3
				th(ii,jj) = sqrt( sum(Wx{ii+1}(times).^2)/N );
			end
			th(4,jj) = .5*sqrt( sum(Wx{ii+1}(times).^2)/N );
		end
		
	%% search for R-peak
		for s = 1%for scale 1
				
				% search for maxima exceeding the thresholds
				IX = [];
				for jj = 1:numBlocks
					times = [1:ceil(N/numBlocks)] + (jj-1)*ceil(N/numBlocks);
					if jj == numBlocks
						times = ((jj-1)*ceil(N/numBlocks)+1):N;
					end
					IX = [IX  (find(abs(Wx{s+1}(times)) > th(s,jj)) + (jj-1)*ceil(N/numBlocks)) ];
				end
					
				% search for zero crossing at scale 2^1
					signs = sign(Wx{s+1}(IX));
					difIX = IX(2:end) - IX(1:end-1);
					difSigns = signs(2:end) - signs(1:end-1);
					jumps = find( (difSigns ~= 0) & (difIX <SEP) ); % find sign changes with a min separation

					% find the closest point to zero
					zeroCrossing = [];
					for jj = 1:numel(jumps)
						[sink tmp(jj)] = min(abs(Wx{2}( IX(jumps(jj)):IX(jumps(jj))+1 ))); 
						zeroCrossing(jj) = IX(jumps(jj))+tmp(jj);
					end
					
					
				% reject isolated and redundant maximum lines
					% implicit before
					
				% apply refractory period
					NZ = numel(zeroCrossing);
					jj =1;
					while jj < NZ

						difZ = zeroCrossing(jj) - zeroCrossing;
						tooClose = false(1,NZ);
						tooClose(abs(difZ) < REFRAC) = 1;
						[sink tmp] = max(abs(x(zeroCrossing).*(tooClose)));
						mask = true(1,NZ);
						mask(tooClose) = 0;
						mask(tmp) = 1;
						zeroCrossing = zeroCrossing(mask);
						NZ = numel(zeroCrossing);
						jj = jj +1;
					end
		end
		Rwave = zeroCrossing;
		
		
		
	%% plot
		if optPlot
				figure;
				hold on;
				title('RMS of the signal');
				plot(Wx{2},'m--');
				line([0 N],[th(1,1) th(1,1)],'Color','r', 'LineStyle','-.');
				line([0 N],[-th(1,1) -th(1,1)],'Color','r', 'LineStyle','-.');
				plot(Wx{3},'c--');
				plot(x);
				plot(Rwave,x(Rwave),'rx','LineWidth',2);
				hold off;
			end

end