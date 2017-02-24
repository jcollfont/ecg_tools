function [alX] = align_signal(X)


	
	[M N] = size(X);
	
	for ii = 1:M;
		for jj = ii:M; 
			corr = xcorr(X(ii,:),X(jj,:));
			[sink ix] = max(abs(corr));
			del(ii,jj) = ix;
		end;
	end;
	
	dif = (N - del(1,:));

	mD = min(dif);
	sN = N - (max(dif) - min(dif));
	ix = [];
	for ii = 1:M;
		ix(ii,:) = [(dif(ii) - mD + 1): (dif(ii) -mD + sN)];
	end;

	alX =zeros(M, sN);
	for ii = 1:M; 
		alX(ii,:) = X(ii,ix(ii,:));
	end;
end