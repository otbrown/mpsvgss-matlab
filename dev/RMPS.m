% RMPS.m
% function which forms the right canonical matrix product state from a state vector
% please note that the function MPS acts as an interface to this one; it's quite okay 
% to use this directly, but worth looking at MPS to see what clean-up it performs
% (also consider just using it...)
% Oliver Thomson Brown
% 03/11/2014
%
% [RETURN]
% matrices	: cell array, A_j_n = matrices{n}(:,:,j)
%
% [INPUT]
% init		: d^N by 1 sparse dbl array, contains initial state vector, where d is the physical dimension and N is the size of the lattice
% HILBY		: int, the physical dimension of a single site
% COMPRESS	: int, declares the SVD compression to be applied to the state, often chi in the literature ***MUST BE 'Inf' FOR EXACT STATE FORMATION***
% L		: int, the number of sites in the system
% SPACE		: the dimension of the total system hilbert space


function [matrices] = RMPS(init, HILBY, COMPRESS, L, SPACE)
	% RETURN ALLOCATION
	matrices = cell(L,1);

	% LAST SITE
	state = transpose( reshape(init,HILBY,SPACE/2) );
	[U, S, V] = Decomp(state, COMPRESS);
	V =  ctranspose(V);

	for sigma = 1:1:HILBY				% site L column vectors
		matrices{L}(:, :, HILBY + 1 - sigma) = V(1:2, sigma);
	end
	
	fprintf('site %d constructed\n', L);		% progress indicator

	chain = U * S;					% all sites to left of L
	route = L-1 : -1 : 2;

	colLim = 2;

	% BULK SITES -- NOT CURRENTLY HILBY > 2 COMPATIBLE!
	for site = route
		[chainRow, chainCol] = size(chain);

		state = zeros(chainRow / 2, 2 * chainCol);
		state( :, 1 : chainCol) = chain( 1 : 2 : end, :);
		state( :, chainCol + 1 : 2 * chainCol) = chain( 2 : 2 : end, :);  		

		rowLim = max(min(sprank( state ), COMPRESS), 2);

		[U, S, V] = Decomp(state, COMPRESS);

		V = ctranspose(V);
		Vcol = size(V,2);

		matrices{site}(:, :, :) = cat(3, V(1:rowLim, Vcol/2 + 1 : Vcol/2 + colLim), V(1:rowLim, 1:colLim));

		chain = U * S;

		colLim = rowLim;

		fprintf('site %d constructed\n', site);		% progress indicator
	end

	% FIRST SITE
	for sigma = 1 : 1 : HILBY
		matrices{1}(:, :, HILBY + 1 - sigma ) = chain(sigma, 1:2);   		
	end

	fprintf('site 1 constructed\n');		% progress indicator
end
