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

function [lmps] = LCan2(mps, route)
    % RETURN ALLOCATION
    lmps = mps;
    HILBY = size(lmps{1}, 3);
    
    for site = route		
        for localState = 1 : 1 : HILBY
            [U, S, V] = svd(lmps{site}(:, :, localState));
            lmps{site}(:, :, localState) = U;

            chain = S * ctranspose(V);
            lmps{site + 1}(:, :, localState) = ...
                chain * lmps{site+1}(:, :, localState);
        end 
    end
end
