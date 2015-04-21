% EffH.m
% function which generates the effective Hamiltonian on a given site, and maps it from a rank-6 tensor to a matrix
% Oliver Thomson Brown
% 2015-03-02
% 
% [RETURN]
% effectiveHamiltonian	: rank-6 tensor reshaped to be a matrix, contains the contraction of the whole tensor network except the update site tensor
% 
% [INPUTS]
% HILBY             : int, dimension of the local state space
% rowMax            : int, the number of rows in the update site tensor
% colMax            : int, the number of columns in the update site tensor
% leftBlock         : rank 3 tensor, contains the contraction through the whole network left of the update site
% mpoTensor			: double array, contains matrix product operator -- mpo{1} of first site, mpo{2} of bulk sites, mpo{3} of last site
% rightBlock		: rank 3 tensor, contains the contraction through the whole network right of the update site  

function [ effectiveHamiltonian ] = EffH(HILBY, rowMax, colMax, leftBlock, mpoTensor, rightBlock)
	% pre-allocate
	dimension = HILBY * rowMax * colMax;
	effectiveHamiltonian = zeros(dimension, dimension);
    
    % permute left and right blocks
    leftBlock = permute(leftBlock, [1,3,2]);
    rightBlock = permute(rightBlock, [3,1,2]);

	% LOOP THE LOOP
    for braState = 0 : 1 : HILBY - 1
        for ketState = 0 : 1 : HILBY - 1
            stateMPO = mpoTensor(braState + 1 : HILBY : end, ketState + 1 : HILBY : end);
            for row = 0 : 1 : rowMax - 1
                for col = 1 : 1 : colMax
                    for conjCol = 0 : 1 : rowMax - 1
                        jRow = braState * colMax * rowMax + conjCol * colMax;
                        jCol = ketState * rowMax * colMax + row * colMax + col;
                        % introduce blocks and contract LWR
                        effectiveHamiltonian(jRow + 1 : jRow + colMax, jCol) = effectiveHamiltonian(jRow + 1 : jRow + colMax, jCol) + ...
                            transpose(leftBlock(conjCol + 1, :, row + 1) * stateMPO * rightBlock(:, :, col));
                        % END TIMES
                    end
                end
            end
        end
    end
	 
end