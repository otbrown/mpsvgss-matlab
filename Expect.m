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
	[rowMax, colMax, HILBY] = size(M);	

	conjM = zeros(colMax, rowMax, HILBY);
	for localState = 1 : 1 : HILBY
		conjM(:, :, localState) = ctranspose( M(:, :, localState) );
	end

	opCount = length(mpo{1}) / HILBY ;
	
	if TARGET == 1			% select correct mpo
		mpodex = 1;
		opRowMax = 0;
		opColMax = opCount - 1;
	elseif TARGET == L
		mpodex = 3;
		opRowMax = opCount - 1;
		opColMax = 0;
	else
		mpodex = 2;
		opRowMax = opCount - 1;
		opColMax = opCount - 1;
	end

	% CALCULATION BEGINS
	expectationValue = 0;
	
	for braState = 1 : 1 : HILBY
		for ketState = 1 : 1 : HILBY
			for row = 1 : 1 : rowMax
				for col = 1 : 1 : colMax
					for conjRow = 1 : 1 : colMax
						for conjCol = 1 : 1 : rowMax
							for opRow = 0 : 1 : opRowMax
								for opCol = 0 : 1 : opColMax
									expectationValue = expectationValue + leftBlock( conjCol, row, opRow + 1) ...
										* mpo{mpodex}( opRow * HILBY + braState, opCol * HILBY + ketState) ...
										* rightBlock( conjRow, col, opCol + 1) * conjM( conjRow, conjCol, braState) ...
										* M( row, col, ketState); 
								end
							end
						end
					end
				end
			end
		end
	end	
end
