% EffH.m
% Oliver Thomson Brown
% 2015-03-02
% 
% [RETURN]
% effectiveHamiltonian	: rank-6 tensor reshaped to be a matrix, contains the contraction of the whole tensor network except the update site tensor
% 
% [INPUTS]
% HILBY			: int, dimension of the local state space
% rowMax		: int, the number of rows in the update site tensor
% colMax		: int, the number of columns in the update site tensor
% leftBlock		: rank 3 tensor, contains the contraction through the whole network left of the update site
% mpo			: cell array, 3 * 1, contains matrix product operator -- mpo{1} of first site, mpo{2} of bulk sites, mpo{3} of last site
% rightBlock		: rank 3 tensor, contains the contraction through the whole network right of the update site  

function [ effectiveHamiltonian ] = EffH(HILBY, rowMax, colMax, leftBlock, mpo, rightBlock);
	% gather data	
	opRowMax = ( size(mpo, 1) / HILBY ) - 1;
	opColMax = ( size(mpo, 2) / HILBY ) - 1;

	% pre-allocate
	dimension = HILBY * rowMax * colMax;
	effectiveHamiltonian = sparse(dimension, dimension);

	% LOOP THE LOOP
	for braState = 0 : 1 : HILBY - 1
		for ketState = 0 : 1 : HILBY - 1
			for row = 0 : 1 : rowMax - 1
				for col = 1 : 1 : colMax
					for conjRow = 1 : 1 : colMax 
						for conjCol = 0 : 1 : rowMax - 1
							% joint indexing
							jRow = braState * colMax * rowMax + conjCol * colMax + conjRow;
							jCol = ketState * rowMax * colMax + row * colMax + col;
							% introduce block indices and contract LWR
							for opRow = 0 : 1 : opRowMax
								for opCol = 0 : 1 : opColMax
									effectiveHamiltonian(jRow, jCol) = effectiveHamiltonian(jRow, jCol) + ...
														leftBlock(conjCol + 1, row + 1, opRow + 1) ...
														* mpo(opRow * HILBY + braState  + 1, opCol * HILBY + ketState + 1) ...
														* rightBlock(conjRow, col, opCol + 1);
									% END TIMES
								end
							end
						end
					end
				end
			end
		end
	end
	 
end
