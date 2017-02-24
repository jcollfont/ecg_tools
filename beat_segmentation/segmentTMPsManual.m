function [stmp t] = segmentTMPsManual(tmp)
%% HELP:
%		segmentTMPs:
%			segments the given TMPs at the user chosen points
%
%		INPUT:
%			- tmp - <M,N>double - matrix containing the tmp series of N
%			samples for the M points on the heart.
%
%		OUTPUT:
%			- stmps - <1,M>cell - each cell contains a cell array with the
%			tmps of a single beat per each point on the heart.
%
%
%

	[M N] = size(tmp);
	
	stmp = cell(1,M);
	t = cell(1,M);
	
	%% PLOT
	fig = figure;
	hold on;
	for m = 1:M
		plot(tmp(m,:));
	end
	hold off;
	
	
	%% CUT
	cuts = 1;
	while (cuts(end) > 0)&(cuts(end)<N)
		[x y] = ginput(1);
		cuts = [cuts floor(x)];
		if (cuts(end) > 0)&(cuts(end)>N)
			cuts(end) = N;
		end
	end
	
	
	nB = numel(cuts);
	%% SEGMENT
	for m = 1:M
		
		stmp{m} = cell(1,nB-1);
		t{m} = cell(1,nB-1);
		
		stmp{m}{1} = tmp(m,1:cuts(2));
		t{m}{1} = 1:cuts(2);
		for b = 2:nB-1
			stmp{m}{b} = tmp(m,cuts(b)+1:cuts(b+1));
			
			t{m}{b} = cuts(b)+1:cuts(b+1);
		end
		
	end
	
	close(fig);
end