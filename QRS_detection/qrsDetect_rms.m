function [QRS] = qrsDetect_rms(x, t)

	
	[L N] = size(x);

	Nw = 50;
	searchT = 100;
	idleT = 100;
	tav = 2;
	tha = 1/3;
	
	rms = sum(x.^2,1);
	
	initMax = max(rms);

	
	[peaks z] = peak_searchPPT(rms, Nw, tha, tav, initMax, searchT, idleT);
	
	QRS = repmat(t(peaks),L,1);
end