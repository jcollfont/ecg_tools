function [tseg qrst beats fullcuts] = segmentHeartBeats_bySections(fullfbPot, nsamplespersec)
% Author: Jaume Coll Font
% Interactively partitions ...etc (TODO: fill this in as we go) 
%
% Usage:

    
    [m n] = size(fullfbPot)
	beats = cell(m,1);
	tseg = cell(m,1);
	qrst = cell(m,1);

	
	numSec = ceil(n/nsamplespersec)
	
	cuts = cell(1,numSec);
	fullcuts = [];
	for sec =1:numSec
		sec
		fbPot = fullfbPot(:,(sec-1)*nsamplespersec + 1: min(sec*nsamplespersec,n) );

		rmsfbPot = sum(fbPot.^2,1);
		
		plot(rmsfbPot);
		
		maxval=max(rmsfbPot(:));
		cuts{sec} = [];
		while true
			[x y button] = ginput(1);
			if(button>1), break; end
			hold on; line(floor(x)*ones(1,2),[0,maxval],'Color','r'), hold off;
			cuts{sec} = [cuts{sec} floor(x)];
			cuts{sec}(end)
			if (cuts{sec}(end) > nsamplespersec)
				break;
			end

		end	
		
		fullcuts = [ fullcuts (cuts{sec}(1:end-1) + (sec-1)*nsamplespersec)];
	end
	fullcuts = [1 fullcuts n];

	for l = 1:m
		count = 1;
		for b = 1:2:numel(fullcuts)-2
			beats{l,count} = fullfbPot(l,fullcuts(b):fullcuts(b+2));
			qrst{l,count} = fullfbPot(l,fullcuts(b):fullcuts(b+1));
			tseg{l,count} = fullfbPot(l,fullcuts(b+1):fullcuts(b+2));
			count = count +1;
		end
	end
	
	
	figure;hold on;
	plot(fullfbPot(1,:),'k--')
	for ii = 1:104
		t = fullcuts(2*ii-1):1:fullcuts(2*ii);
		plot(t,qrst{1,ii})
		t = fullcuts(2*ii):1:fullcuts(2*ii+1);
		plot(t,tseg{1,ii},'r')
	end
	plot(fullcuts,0,'go');
	hold off;

end