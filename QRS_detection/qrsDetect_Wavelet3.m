function [QRS] = qrsDetect_Wavelet3(x, t, optPlot)
%%	HELP:
%		[QRS_peak freq] = qrsDetect_Wavelet3(x,t)
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
		QRS = cell(1,L);
		
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
					Wx{ii+1} = tmp;
				end

		%% define the thresholds
			th = zeros(1,4);
			for ii =1:3
				th(ii) = sqrt( sum(Wx{ii+1}.^2)/N );
			end
				th(4) = .5*sqrt( sum(Wx{ii+1}.^2)/N );
			
		%% search for QRS peak
			for s = 1%for scale 1
				
				% search for maxima exceeding the thresholds
					IX = find(abs(Wx{s+1}) > th(s));
					
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
						[sink tmp] = max(abs(x(l,zeroCrossing).*(tooClose)));
						mask = true(1,NZ);
						mask(tooClose) = 0;
						mask(tmp) = 1;
						zeroCrossing = zeroCrossing(mask);
						NZ = numel(zeroCrossing);
						jj = jj +1;
					end

				% search back with lower thresholds if no peak found in
				% interval.
					% "pa luego"
		
			end
		
			QRS_peak = zeroCrossing;
		
		%% search for QRS individual waves
		
			% define search window
				sw = 10;	%%%%%%%%%%%%%%%%%%%%%%%%%HARDCODED, MIGHT VARY DEPENDING ON SAMPLING RATE AND BEAT RATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% define thresholds
				Gpre = 0.06*max(Wx{3});
				Gpost = 0.09*max(Wx{3});
			
			% for each qrs peak 
			Nqrs = numel(QRS_peak);
			QRS{l} = zeros(Nqrs,3);
			for jj = 1:Nqrs
				
				% find the moduli around qrs
					start = max((QRS_peak(jj)-floor(SEP/2)),1);
					finish = min((QRS_peak(jj)+floor(SEP/2)),N);

					[vpre npre] = max(abs(   Wx{3}(start:QRS_peak(jj))  ));
					[vpost npost] = max(abs( Wx{3}((QRS_peak(jj)):finish) ));
					vpre = vpre + Wx{3}(QRS_peak(jj));
					vpost = vpost + Wx{3}(QRS_peak(jj));
					npre = start + npre -1;
					npost = QRS_peak(jj) + npost -1;
				
				% find other peaks around the moduli
					start = max( npre - sw, 1 );
					finish = min( npost + sw, N );
					
					IXpre = find( abs(Wx{3}(start:npre)) > Gpre ) + start -1;
					IXpost = find( abs(Wx{3}(npost:finish)) > Gpost ) + npost -1;
					
					difIXpre = sign(Wx{3}(IXpre(2:end))) + sign(Wx{3}(IXpre(1:end-1)));
					difIXpost = sign(Wx{3}(IXpost(2:end))) + sign(Wx{3}(IXpost(1:end-1)));
					

					jumpspre = find(difIXpre == 0);
					jumpspost =  find(difIXpost == 0);
			
					
					
				% find zero crossing between peaks that exceed the
				% threshold
					zeroPre = [];
					if (numel(IXpre) ~= 0) & ( numel(jumpspre) == 1 )
						[sink zeroPre] = min(abs( Wx{2}(IXpre(jumpspre(1):(jumpspre(1)+1))) ));
						zeroPre = IXpre(jumpspre) + zeroPre -1;
					elseif (numel(IXpre) ~= 0) & ( numel(jumpspre) > 1 )
						[sink zeroPre(1)] = min(abs( Wx{2}(IXpre(jumpspre(1):(jumpspre(1)+1))) ));
						zeroPre(1) = IXpre(jumpspre(1)) + zeroPre(1) -1;
						[sink zeroPre(2)] = min(abs( Wx{2}(IXpre(jumpspre(2):(jumpspre(2)+1))) ));
						zeroPre(2) = IXpre(IXpre(jumpspre(2))) + zeroPre -1;
					end
					
					zeroPost = [];
					if (numel(IXpost) ~= 0) & ( numel(jumpspost) == 1 )
						[sink zeroPost] = min(abs( Wx{2}(npost:IXpost(jumpspost(1))) ));
						zeroPost = IXpost(jumpspost) + zeroPost -1;
					elseif (numel(IXpost) ~= 0) & ( numel(jumpspost) > 1 )
						[sink zeroPost(1)] = min(abs( Wx{2}(IXpost(jumpspost(1):(jumpspost(1)+1))) ));
						zeroPost(1) = IXpost(jumpspost(1)) + zeroPost(1) -1;
						[sink zeroPost(2)] = min(abs( Wx{2}(IXpost(jumpspost(2):(jumpspost(2)+1))) ));
						zeroPost(2) = IXpost(jumpspost(2)) + zeroPost(2) -1;
					end
					
				% find structure QRS, QR, RS
					type = 1*numel(zeroPre) + 10*numel(zeroPost);
					switch type
						case 0
							QRS{l}(jj,2) = QRS_peak(jj);
							
						case 1
							QRS{l}(jj,1:2) = [zeroPre(1) QRS_peak(jj)];
							
						case 2
							QRS{l}(jj,1:3) = [zeroPre(1) zeroPre(2) QRS_peak(jj)];
							
						case 10
							QRS{l}(jj,2:3) = [ QRS_peak(jj) zeroPost(1) ];
						
						case 20 
							QRS{l}(jj,1:3) = [ QRS_peak(jj) zeroPost(1) zeroPost(2)];
							
						otherwise
							QRS{l}(jj,1:3) = [ zeroPre(1) QRS_peak(jj) zeroPost(1) ];
								
					end
					
			end

		%% plot
			if optPlot
				figure;
				hold on;
				title(sprintf('Lead: %d',l));
				plot(Wx{2},'m--');
				line([0 N],[th(1) th(1)],'Color','r', 'LineStyle','-.');
				line([0 N],[-th(1) -th(1)],'Color','r', 'LineStyle','-.');
				plot(Wx{3},'c--');
				line([0 N],[Gpre Gpre],'Color','c', 'LineStyle','-.');
				line([0 N],[-Gpre -Gpre],'Color','c', 'LineStyle','-.');
				line([0 N],[Gpost Gpost],'Color','g', 'LineStyle','-.');
				line([0 N],[-Gpost -Gpost],'Color','g', 'LineStyle','-.');
				plot(x(l,:));
				plot(QRS{l}(QRS{l}(:,1)~=0,1),x(l,QRS{l}(QRS{l}(:,1)~=0,1)),'go','LineWidth',2);
				plot(QRS{l}(QRS{l}(:,3)~=0,3),x(l,QRS{l}(QRS{l}(:,3)~=0,3)),'ko','LineWidth',2);
				plot(QRS{l}(QRS{l}(:,2)~=0,2),x(l,QRS{l}(QRS{l}(:,2)~=0,2)),'rx','LineWidth',2);
				hold off;
			end
	end
	
end