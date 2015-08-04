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
	
    HILBY = size(lmps{1}, 3);    

	for site = route		
        [rowSz, colSz, ~] = size(lmps{site});

		M = reshape(permute(lmps{site}, [1, 3, 2]), [HILBY * rowSz, colSz]);

        [Q, R] = qr(M, 0);

        lmps{site} = permute(reshape(Q, [rowSz, HILBY, colSz]), [1, 3, 2]);

        for localState = 1 : 1 : HILBY
            lmps{site + 1}(:, :, localState) = R * lmps{site + 1}(:, :, localState);
        end 
	end
end
