% Expect.m
% function to calculate the expectation value of an MPO over an MPS
% Oliver Thomson Brown
% 15-01-29
% DOCSTRING!

function [ expectationValue ] = Expect(mps, mpo, leftBlock, rightBlock, TARGET)

	% TASHA YAR
	L = size(mps, 1);

	M = mps{TARGET};
	[rowMax, colMax, HILBY] = size(M);

	conjM = zeros(colMax, rowMax, HILBY);
	for localState = 1 : 1 : HILBY
		conjM(:, :, localState) = ctranspose( M(:, :, localState) );
	end

	opCount = ( length(mpo{1} / HILBY ) - 1;
	
	if TARGET == 1
		mpodex = 1;
		opRowMax = 0;
		opColMax = opCount;
	elseif TARGET == L
		mpodex = 3;
		opRowMax = opCount;
		opColMax = 0;
	else
		mpodex = 2;
		opRowMax = opCount;
		opColMax = opCount;
	end

	% CALCULATION BEGINS
	expectationValue = 0;
	
	for braState = 1 : 1 : HILBY
		for ketState = 1 : 1 : HILBY
			for row = 1 : 1 : rowMax
				for col = 1 : 1 : colMax
					for conjRow = 1 : 1 : colMax
						for conjCol = 1 : 1 : rowMax
							for opRow = 1 : 1 : opRowMax
								for opCol = 1 : 1 : opColMax
									expectationValue = expectationValue + leftBlock(conjRow, opRow + 1, row) ...
										* mpo{mpodex}(opRow * HILBY + braState, opCol * HILBY + ketState) ...
										* rightBlock(conjCol, opCol + 1, col) * conjM(conjRow, conjCol, braState) ...
										* M(row, col, ketState);
								end
							end
						end
					end
				end
			end
		end
	end

end
