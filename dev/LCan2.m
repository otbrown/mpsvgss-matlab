% LCan.m
% function which left-canonises an already formed matrix product state
% please note that the function Can acts as an interface to this one
% in particular Can ensures that the supplied route is valid 
% Oliver Thomson Brown
% 10/11/2014
%
% [RETURN]
% lmps		: cell array containing left-canonised MPS, A_state_site = lmps{site}(:,:,state)
%
% [INPUTS]
% mps		: cell array containing a matrix product state, A_state_site = mps{site}(:,:,state)
% route		: vector containing sites which are to be canonised, if whole chain is to be canonised 
%		  route should be 1 : 1 : L-1 where L is the length of the chain

function [lmps] = LCan(mps, route)
	% RETURN ALLOCATION
	lmps = mps;
	
    [rowMax, colMax, HILBY] = size(lmps{route(1)});
    
	for site = route		
		perm = permute(lmps{site}, [1, 3, 2]);
		M = reshape(perm, [HILBY * rowMax, colMax]);
		[U, S, V] = svd(M);

		U = U(1 : (HILBY * rowMax), 1 : colMax);
		reshU = reshape(U, [rowMax, HILBY, colMax]);
		lmps{site} = permute(reshU, [1, 3, 2]);

		chain = S * ctranspose(V);

		rowMax = colMax;
		colMax = size(lmps{site + 1}(:,:,1), 2);
		
		for localState = 1 : 1 : HILBY
			N = chain * lmps{site+1}(:, :, localState);
			lmps{site + 1}(:, :, localState) = N(1 : rowMax, 1 : colMax);
		end 
	end
end
