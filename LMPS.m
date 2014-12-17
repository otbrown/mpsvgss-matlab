% LMPS.m
% function which forms the left canonical matrix product state from a state vector
% please note that the function MPS acts as an interface to this one; it's quite okay 
% to use this directly, but worth looking at MPS to see what clean-up it performs
% (also consider just using it...)
% Oliver Thomson Brown
% 16/10/2014
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

function [matrices] = LMPS(init, HILBY, COMPRESS, L, SPACE)	
	
	threshold = 1E-14;				% mid-calculation values below this threshold will be zeroed -- improves stability

	% RETURN ALLOCATION
	matrices = cell(L,1);

	% FIRST SITE
	state = transpose( reshape(init, SPACE/HILBY, HILBY) );
	[U, S, V] = Decomp(state, COMPRESS);		% Decomp filters exact states to svd(full(state)) which is faster than svds for all the singular values of a matrix

	for sigma = 1 : HILBY
		matrices{1}(:, :, (HILBY - sigma + 1) ) = U(sigma, :);   		% site 1 vectors
	end

	fprintf('site 1 constructed\n');

	chain = S * ctranspose(V);			% all sites to right of site 1
	route = 2 : 1 : L-1;

	rowLim = 2;

	% BULK SITES -- NOT CURRENTLY HILBY > 2 COMPATIBLE!
	for site = route

		chainCol = size(chain,2);

		state = cat(1,chain(:, 1:chainCol/2), chain(:, chainCol/2+1 : chainCol));
		
		colLim = max(min(sprank( state ), COMPRESS), 2);

		[U, S, V] = Decomp(state, COMPRESS);

		U(abs(U) < threshold) = 0;
		S(abs(S) < threshold) = 0;
		V(abs(V) < threshold) = 0;

		Urow = size(U,1);

		matrices{site}(:, :, :) = cat(3, U(Urow/2 + 1 : Urow/2 + rowLim, 1:colLim), U(1:rowLim, 1:colLim));

		chain = S * ctranspose(V);

		fprintf('site %d constructed\n', site);

		rowLim = colLim;
	end

	% LAST SITE
	for sigma = 1 : HILBY
		matrices{L}(:, :, (HILBY - sigma + 1) ) = chain(1:2, sigma);   		
	end

	fprintf('site %d constructed\n', L);		% progress indicator

end
