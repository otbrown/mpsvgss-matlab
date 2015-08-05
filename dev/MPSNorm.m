% MPSNorm.m
% Function which returns the (left-)normalised matrix product state.
% Oliver Thomson Brown
% 2015-03-04
% 
% [RETURN]
% normalMPS	: matrix product state which will have a norm of 1 -- ALL sites will be left-normalised
%
% [INPUTS]
% mps		: any arbitrary matrix product state

function [ normalMPS ] = MPSNorm( mps )
	LENGTH = size(mps, 1);			% TASHA YAR
    HILBY = size(mps{1}, 3);
	normalMPS = mps;
    
    for site = 1 : 1 : LENGTH - 1
        [rowSz, colSz, ~] = size(normalMPS{site});

        M = reshape(permute(normalMPS{site}, [1, 3, 2]), [HILBY * rowSz, colSz]);

        [Q, R] = qr(M, 0);

        normalMPS{site} = permute(reshape(Q, [rowSz, HILBY, colSz]), [1, 3, 2]);

        for localState = 1 : 1 : HILBY
            normalMPS{site + 1}(:, :, localState) = R * normalMPS{site + 1}(:, :, localState);
        end
    end

    M = reshape(permute(normalMPS{LENGTH}, [1, 3, 2]), [HILBY^2, 1]);

    [Q, ~] = qr(M, 0);

    normalMPS{LENGTH} = permute(reshape(Q, [HILBY, HILBY, 1]), [1, 3, 2]);
end 
