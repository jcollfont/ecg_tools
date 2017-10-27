function [allPots] = interpolateMissingData( geom, leadlinks, measuredPots)

	[Lapl] = computeLaplacianMatrix(geom,1);
	
	N = size(Lapl,1);
	Nk = numel(find(leadlinks));

	leadlinksBool = false(N,1);
	leadlinksBool(leadlinks) = true;
	
	Lu = Lapl(~leadlinksBool,~leadlinksBool);
	Lk = Lapl(leadlinksBool,leadlinksBool);
	Lc = Lapl(~leadlinksBool,leadlinksBool);

% 	interpPots = (eye(size(Lu)) - Lu)\ Lc * measuredPots;
	interpPots = -(Lu'*Lu + Lc*Lc') \ ( Lc*Lk + Lu*Lc ) * measuredPots;

% 	L1 = [ Lk ; Lc ];
% 	L2 =  [ Lc' ; Lu ];
% 	interpPots = -(L2'*L2) \ L2'*L1 * measuredPots;
	
	allPots = zeros(N,size(measuredPots,2));
	
	allPots(leadlinksBool,:) = measuredPots;
	allPots(~leadlinksBool,:) = interpPots;
	
end