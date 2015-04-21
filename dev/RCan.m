% RCan.m
% function which right-canonises an already formed matrix product state
% please note that the function Can acts as an interface to this one
% in particular Can ensures that the supplied route is valid 
% Oliver Thomson Brown
% 13/11/2014
%
% [RETURN]
% rmps		: cell array containing right-canonised MPS, A_state_site = lmps{site}(:,:,state)
%
% [INPUTS]
% mps		: cell array containing a matrix product state, A_state_site = mps{site}(:,:,state)
% route		: vector containing sites which are to be canonised, if whole chain is to be canonised 
%		  route should be L : -1 : 2 where L is the length of the chain

function [rmps] = RCan(mps, route)
	% RETURN ALLOCATION
	rmps = mps;
    
    [rowMax, colMax, HILBY] = size(rmps{route(1)});
	
    for site = route		
        M = reshape(rmps{site}, [rowMax, colMax * HILBY]);
		[U, S, V] = svd(M);

		V = ctranspose(V);
        V = V(1 : rowMax, 1 : (colMax * HILBY));

        rmps{site} = reshape(V, [rowMax, colMax, HILBY]);

		chain = U * S;

        colMax = rowMax;
        rowMax = size(rmps{site - 1}(:,:,1), 1); 

        for localState = 1 : 1 : HILBY
            N = rmps{site - 1}(:, :, localState) * chain;
            rmps{site - 1}(:, :, localState) = N(1 : rowMax, 1 : colMax);
        end
    end
end
