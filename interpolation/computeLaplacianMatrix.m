function [Lapl] = computeLaplacianMatrix(geom,pathLength)

	if ~exist('pathLength')
		pathLength = 2;
	end
	
	[AdjMtrx] = computeAdjacencyMatrix(geom, pathLength);
	
	Lapl = diag(sum(AdjMtrx,2)) - AdjMtrx;
	

end