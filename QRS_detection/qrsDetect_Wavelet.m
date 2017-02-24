function [QRS_peak freq Wx] = qrsDetect_Wavelet(x, t, optPlot)
%%	HELP:
%		[QRS_peak freq] = qrsDetect_Wavelet(x,t)
%			This algorithm implements a qrs detection algorithm based on
%			the wavelet transform.
%			The reference paper used is:
%			Juan Pablo Martínez et al., "A Wavelet-Based ECG Delineator:
%			Evaluation on Standard Databases", IEEE transactions on
%			biomedical engineering, April 2004.
%
%		INPUT:
%			- x - <L,N>double - input ECG.
%			- t - <1,N>double - sampling times of the ECG.
%			- optPlot - boolean - true if the user wants to plot the
%			results.
%
%		OUTPUT:
%			- QRS_peak - <1,JJ>cell - each cell contains the 0 crossing
%				points measured in each of the scales. These are ordered from s
%				= 2^1 to 2^JJ.
%			- freq - <L,JJ>double - frequency of the signal. This has been
%				measured as the inverse of the average period between QRS'.
%					
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
		SEP = 5; %max separation between jumps
		JJ = 10; % num of sampled scales
		[L N] = size(x);
		
		if ~exist('optPlot')
			optPlot = false;
		end
			

	%% Apply wavelet transform
		h = [1 3 3 1 0]/8;
		g = [0 1 -1 0 0]*2;

	%% for all leads
	for l= 1:L
			
			Wx = {x(l,:)};Sx={x(l,:)};for jj = 1:JJ; Wx{jj+1} = filter(g,1,Sx{jj}(1:2:end));Sx{jj+1} = filter(h,1,Sx{jj}(1:2:end));end;

		%% define the thresholds
			th = zeros(1,4);
			for ii =1:3
				th(ii) = sqrt( sum(Wx{ii+1}.^2)/N );
			end
				th(4) = .5*sqrt( sum(Wx{ii+1}.^2)/N );

		%% search for QRS peak
			for ii = 1%1:4

				% find the waves
				IX{ii} = find(abs(Wx{ii+1}) > th(ii));

				% find sign changes in the wavelet
				signs = sign(Wx{ii+1}(IX{ii}));
				difIX = IX{ii}(2:end) - IX{ii}(1:end-1);
				difSigns = signs(2:end) - signs(1:end-1);
				jumps = find( (difSigns ~= 0) & (difIX <SEP) );

				% find the 0 crossing in between
				for jj = 1:numel(jumps)
					[sink tmp(jj)] = min(abs(Wx{ii+1}( IX{ii}(jumps(jj)):IX{ii}(jumps(jj))+1 ))); 
					zeroCrossing{ii}(jj) = IX{ii}(jumps(jj))+tmp(jj) - 2*ii; 
				end

				% correct the subsampling
				mask = zeros(1,numel(Wx{ii+1}));
				for jj = 1:ii
					NZ = numel(mask);
					mask = zeros(1,NZ);
					mask(zeroCrossing{ii}) = 1;
					mask = reshape( [mask; zeros(1,NZ)],1,2*NZ );
				end

				% set peaks
				QRS_peak{ii}(l,:) = find(mask);


				% measure beat frequency
				freq(l,ii) = mean( 1./(t(QRS_peak{ii}(2:end)) - t(QRS_peak{ii}(1:end-1))) );

			end
			
			if optPlot
				figure;plot(x(l,:));hold on;plot(QRS_peak{ii}(l,:),x(l,QRS_peak{ii}(l,:)),'ro');hold off;
			end
	end
	
end
