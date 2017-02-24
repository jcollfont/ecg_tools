function [QRS_peak freq Wx QRS] = qrsDetect_Wavelet2(x, t, optPlot)
%%	HELP:
%		[QRS_peak freq] = qrsDetect_Wavelet2(x,t)
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
%			- QRS - <Nqrs, 3>double - indicates the index position of Q,R and S
%			waves in this order. If one of the side waves does not exist,
%			it is marked as 0.
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
		SEP = 20; %max separation between jumps
		JJ = 10; % num of sampled scales
		REFRAC = 100; % refractory period
		[L N] = size(x);
		
		if ~exist('optPlot')
			optPlot = false;
		end
			

	%% Apply wavelet transform
		h = [1 3 3 1 0]/8;
		g = [0 1 -1 0 0]*2;

	%% for all leads
	for l= 1:L
			
		%% calculate wavelet transforms
			Wx = {x(l,:)};Sx={x(l,:)};for jj = 1:JJ; Wx{jj+1} = filter(g,1,Sx{jj}(1:2:end));Sx{jj+1} = filter(h,1,Sx{jj}(1:2:end));end;
			
				% correct subsampling in WT
				for ii = 1:JJ
					tmp = Wx{ii+1};
					for jj = 1:ii
						NZ = numel(tmp);
						tmp = reshape(repmat(tmp,2,1),1,2*NZ);
						tmp = [tmp(4:end) zeros(1,4-1)];
					end
					Wx{ii+1} = tmp;%[tmp(4*(ii):N) zeros(1,4*(ii)-1)];
				end

		%% define the thresholds
			th = zeros(1,4);
			for ii =1:3
				th(ii) = sqrt( sum(Wx{ii+1}.^2)/N );
			end
				th(4) = .5*sqrt( sum(Wx{ii+1}.^2)/N );
			
		%% search for QRS peak
			% find the waves
				IX = find(abs(Wx{2}) > th(1));
		
				% find sign changes in the wavelet
				signs = sign(Wx{2}(IX));
				difIX = IX(2:end) - IX(1:end-1);
				difSigns = signs(2:end) - signs(1:end-1);
				jumps = find( (difSigns < 0) & (difIX <SEP) );
								
% 				plot(Wx{2});hold on; line([1 N], [th(1) th(1)],'Color','m'); line([1 N], [-th(1) -th(1)],'Color','m'); plot(IX(jumps),Wx{2}(IX(jumps)),'ro');hold off;pause;
				
				% find the 0 crossing in between
				zeroCrossing = [];
				for jj = 1:numel(jumps)
					[sink tmp(jj)] = min(abs(Wx{2}( IX(jumps(jj)):IX(jumps(jj))+1 ))); 
					zeroCrossing(jj) = IX(jumps(jj))+tmp(jj);
				end
				
% 				plot(Wx{2});hold on; line([1 2500], [th(1) th(1)],'Color','m'); line([1 2500], [-th(1) -th(1)],'Color','m'); plot(zeroCrossing,Wx{2}(zeroCrossing),'ro'); plot(IX(jumps),Wx{2}(IX(jumps)),'kx','LineWidth',2);hold off;pause;
	
				NZ = numel(zeroCrossing);
				if l>1
				if NZ < numel(QRS_peak(l-1,:))
					jumps = find( (difSigns > 0) & (difIX <SEP) );
					zeroCrossing = [];
					for jj = 1:numel(jumps)
						[sink tmp(jj)] = min(abs(Wx{2}( IX(jumps(jj)):IX(jumps(jj))+1 ))); 
						zeroCrossing(jj) = IX(jumps(jj))+tmp(jj);
					end
				end
				end
					
				NZ = numel(zeroCrossing);
				jj =1;
				while jj < NZ
					
					difZ = zeroCrossing(jj) - zeroCrossing;
					tooClose = false(1,NZ);
					tooClose(abs(difZ) < REFRAC) = 1;
					[sink tmp] = max(abs(x(l,zeroCrossing).*(tooClose)));
					mask = true(1,NZ);
					mask(tooClose) = 0;
					mask(tmp) = 1;
					zeroCrossing = zeroCrossing(mask);
					NZ = numel(zeroCrossing);
					jj = jj +1;
				end
				
				l
				NZ
				QRS_peak(l,:) = zeroCrossing;

				freq(l) = mean(mean( 1./(t(QRS_peak(l,2:end)) - t(QRS_peak(l,1:end-1))) ));
	
			
			
		%% FIND INSTANCES Q,R,S
		
			% define search window
			sw = 70;	%%%%%%%%%%%%%%%%%%%%%%%%%HARDCODED, MIGHT VARY DEPENDING ON SAMPLING RATE AND BEAT RATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% define thresholds
			Gpre = .06*max(Wx{3});
			Gpost = 0.09*max(Wx{3});
			
			Nqrs = numel(QRS_peak(l,:));
			for jj = 1:Nqrs
				
				% search for moduli around qrs
				start = max((QRS_peak(l,jj)-sw),1);
				finish = min((QRS_peak(l,jj)+sw),N);
				
				[vpre npre] = max(abs(   Wx{3}(start:QRS_peak(l,jj))   -  (sign(Wx{3}(QRS_peak(l,jj))) == -1)*Wx{3}(QRS_peak(l,jj))   ));
				[vpost npost] = max(abs( Wx{3}((QRS_peak(l,jj)):finish) - (sign(Wx{3}(QRS_peak(l,jj))) == 1)*Wx{3}(QRS_peak(l,jj)) ));
				vpre = vpre + Wx{3}(QRS_peak(l,jj));
				vpost = vpost + Wx{3}(QRS_peak(l,jj));
				npre = start + npre -1;
				npost = QRS_peak(l,jj) + npost -1;
			
				% define type of wave (qrs, rsr', qr, rs, r or qs)
				moduli(l,jj,:) = [ npre , QRS_peak(l,jj) , npost ];
							
% 				plot(Wx{2});hold on;plot(Wx{3},'g--');plot(start:finish,Wx{2}(start:finish),'r'); plot(moduli(l,jj,1),Wx{3}(moduli(l,jj,1)),'go');plot(moduli(l,jj,2),Wx{3}(moduli(l,jj,2)),'kx');plot(moduli(l,jj,3),Wx{3}(moduli(l,jj,3)),'mo');hold off;
				
				% find 0 corssing
				start = max(moduli(l,jj,1)-sw,1);
				finish = min(moduli(l,jj,3)+sw,N);
				s = sign( Wx{2}( start:moduli(l,jj,1) ) );
				ind = find((s(2:end) - s(1:end-1)) ~= 0); % maybe apply refractory period (somehow done) and also eliminate sudden peaks
				if numel(ind ~= 0)
					p1 = ind(end);
				else
					[sink p1] = min(abs( Wx{2}(start:moduli(l,jj,1)) ));
				end
				
				s = sign( Wx{2}( moduli(l,jj,3):finish ));
				ind = find((s(2:end) - s(1:end-1)) ~= 0);
				if numel(ind ~= 0)
					p2 = ind(1);
				else
					[sink p2] = min(abs( Wx{2}(moduli(l,jj,1)):finish ));
				end
				
				p1 = p1 + start -1;
				p2 = p2 + moduli(l,jj,3) -1;
				
				QRS(l,jj,:) = [p1*(vpre > Gpre) moduli(l,jj,2) p2*(vpost > Gpost)];
				
			end			
			
			
			if optPlot
				figure;
				hold on;
				title(sprintf('Lead: %d',l));
				plot(Wx{2},'m--');
				line([0 N],[th(1) th(1)],'Color','r', 'LineStyle','-.');
				line([0 N],[-th(1) -th(1)],'Color','r', 'LineStyle','-.');
				plot(Wx{3},'c--');
				line([0 N],[th(1) th(1)],'Color','c', 'LineStyle','-.');
				line([0 N],[-th(1) -th(1)],'Color','c', 'LineStyle','-.');
				plot(x(l,:));
				plot(QRS(l,QRS(l,:,1)~=0,1),x(l,QRS(l,QRS(l,:,1)~=0,1)),'go','LineWidth',2);
				plot(QRS(l,QRS(l,:,2)~=0,2),x(l,QRS(l,QRS(l,:,2)~=0,2)),'ro','LineWidth',2);
				plot(QRS(l,QRS(l,:,3)~=0,3),x(l,QRS(l,QRS(l,:,3)~=0,3)),'ko','LineWidth',2);
				hold off;
			end
	end
	
end
