function [tseg qrst beats cuts] = segmentHeartBeats(fbPot)
    
    [m n] = size(fbPot);
    beats = {};
	tseg = cell(1,m);
	qrst = cell(1,m);
    
    rmsfbPot = sum(fbPot.^2,1);
    
    plot(rmsfbPot);
	cuts = 1;
	while (cuts(end) > 0)&(cuts(end)<n)
		[x y] = ginput(1);
		cuts = [cuts floor(x)];
		if (cuts(end) > 0)&(cuts(end)>n)
			cuts(end) = n;
		end
	end
       
	
    for l = 1:m
        for b = 2:numel(cuts)
            beats{l}{b-1} = fbPot(l,cuts(b-1):cuts(b));
			
			if mod(b,2) == 0
				qrst{l}{ceil(b/2)} = fbPot(l,cuts(b-1):cuts(b));
			else
				tseg{l}{floor(b/2)} = fbPot(l,cuts(b-1):cuts(b));
			end
		end
		
		
	end
	

	
	
end