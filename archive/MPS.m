% MPS.m
% function which returns the matrix product form of the input state
% exact or compressed, dependent on input parameter COMPRESS
% interfaces LMPS.m for left-canonical MPS and RMPS.m for right-canonical MPS 
% Oliver Thomson Brown
% 16/10/2014
%
% [RETURN]
% matrices	: cell array, A_j_n = matrices{n}(:,:,j), where j indexes state, n indexes site
%
% [INPUT]
% init		: d^N by 1 sparse dbl array, contains initial state vector, where d is the physical dimension and N is the size of the lattice
% HILBY		: int, the physical dimension of a single site
% COMPRESS	: int, declares the SVD compression to be applied to the state, 0 flags an exact MPS, often chi in the literature
% DIRECTION	: char, declares left or right canonical MPS formation (should be 'L' or 'R')

function [matrices] = MPS(init, HILBY, COMPRESS, DIRECTION)

	% INPUT TIDY-UP

	if DIRECTION ~= 'L'
		if DIRECTION ~= 'R'
			disp('Direction invalid! Please enter L for left-canonical matrix product state formation, or R for right-canonical matrix product state formation.');
			return;
		end
	end

	SPACE = size(init, 1);					% determine full state space, SPACE, and length of chain, L
	L = log(SPACE) / log(HILBY);

	init = init / sqrt(ctranspose(init) * init);		% normalise the state vector

	if COMPRESS == 0					% COMPRESS == 0 flags exact state decomposition
		COMPRESS = Inf;
	end

	% DECOMPOSITION INTERFACE

	if DIRECTION == 'L'
		matrices = LMPS(init, HILBY, COMPRESS, L, SPACE);
	elseif DIRECTION == 'R'
		matrices = RMPS(init, HILBY, COMPRESS, L, SPACE);
	end

end
