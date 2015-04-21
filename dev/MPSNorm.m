% MPSNorm.m
% Function which returns the (left-)normalised matrix product state.
% Oliver Thomson Brown
% 2015-03-04
% 
% [RETURN]
% normalMPS	: matrix product state which will have a norm of 1 -- ALL sites will be left-normalised
%
% [INPUTS]
% mps		: any arbitrary matrix product state with a local Hilbert space dimension of 2 (this should be changed one day maybe)

function [ normalMPS ] = MPSNorm( mps )
	L = size( mps, 1 );			% TASHA YAR
	normalMPS = mps;
    
    [rowMax, colMax, HILBY] = size(normalMPS{1});

	for site = 1 : 1 : L - 1	
		perm = permute(normalMPS{site}, [1, 3, 2]);
		M = reshape(perm, [HILBY * rowMax, colMax]);
		[U, S, V] = svd(M);

		U = U(1 : (HILBY * rowMax), 1 : colMax);
		reshU = reshape(U, [rowMax, HILBY, colMax]);
		normalMPS{site} = permute(reshU, [1, 3, 2]);

		chain = S * ctranspose(V);

		rowMax = colMax;
		colMax = size(normalMPS{site + 1}(:,:,1), 2);
		
		for localState = 1 : 1 : HILBY
			N = chain * normalMPS{site+1}(:, :, localState);
			normalMPS{site + 1}(:, :, localState) = N(1 : rowMax, 1 : colMax);
		end 
	end

        [rowMax, colMax, HILBY] = size(normalMPS{L});
		
		perm = permute(normalMPS{L}, [1, 3, 2]);
		M = reshape(perm, [HILBY * rowMax, colMax]);
		[U, S, V] = svd(M);

		U = U(1 : (HILBY * rowMax), 1 : colMax);
		reshU = reshape(U, [rowMax, HILBY, colMax]);
		normalMPS{L} = permute(reshU, [1, 3, 2]);
end 
