% Expect.m
% function to calculate the expectation value of a mixed canonical MPO over an MPS
% currently only allows for the bra and ket mps to be the same, but this will be changed
% Oliver Thomson Brown
% 15-01-29
% 
% [RETURN]
% expectationValue	: scalar double, contains the full contraction through a tensor network including a matrix product operator
%
% [INPUTS]
% mps			: cell array, contains the matrix product state
% mpo			: 3 * 1 cell array, contains the matrix product operator -- mpo{1} first site mpo, mpo{2} bulk site mpo, mpo{3} last site mpo
% leftBlock		: contains the tensor network contraction from the first site up to the target site
% rightBlock		: contains the tensor network contraction from the last site up to the target site
% TARGET 		: int, the target site -- the mps should be left-normalised from the first site to this one, and right-normalised from the last site to this one

function [ expectationValue ] = Expect(mps, mpo, leftBlock, rightBlock, TARGET)

	% TASHA YAR
	L = size(mps, 1);

	M = mps{TARGET};
    conjM = conj(permute(M, [2, 1, 3]));
	[rowMax, colMax, HILBY] = size(M);
    
    if TARGET == 1			% select correct mpo
        mpodex = 1;
    elseif TARGET == L
        mpodex = 3;
    else
        mpodex = 2;
    end

    mpoTensor = mpo{mpodex};
    
    % permutations
    lBlock = permute(leftBlock, [1, 3, 2]);
    rBlock = permute(rightBlock, [3, 1, 2]);
    
	% CALCULATION BEGINS
	expectationValue = 0;
	
    for braState = 1 : 1 : HILBY
        for ketState = 1 : 1 : HILBY
            stateMPO = mpoTensor(braState : HILBY : end, ketState : HILBY : end);
            for row = 1 : 1 : rowMax
                for col = 1 : 1 : colMax
                    expectationValue = expectationValue + sum( diag( ...
                        lBlock( :, :, row) * stateMPO * rBlock( :, :, col) ... 
                        * conjM( :, :, braState) )) * M( row, col, ketState);
                end
            end
        end
    end
end
