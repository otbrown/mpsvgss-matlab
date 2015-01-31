% RBlock.m
% Oliver Thomson Brown
% 2015-01-15
%
% [RETURN]
% rightBlock	: 3-d double array which contains contraction of sites right of the target site 
%
% [ARGUMENTS]
% mps		: L by 1 cell array containing matrix product state for the system, indexed as mps{site}(:,:,local_state)
% mpo		: 3 by 1 cell array containing matrix product operator, mpo{1} for first site, mpo{2} for bulk sites, mpo{3} for last site
% TARGET	: double, integer value -- target site for update during variational search
%
% *** IMPORTANT ***
% For ease of referencing operators are ordered with <0|operator|0> as the top left element!

function [ rightBlock ] = RBlock(mps, mpo, TARGET)

	% PULL DATA LIKE TASHA YAR
	L = size(mps, 1);
	HILBY = size(mps{L}, 3);
	OPCOUNT = size(mpo{3}, 1) / HILBY;

	inner = ones(2, OPCOUNT, 2);

	for site = L : -1 : TARGET + 1
		if site == L		% sets matrix product operator choice and sizes
			mpodex = 3;
			opRowMax = OPCOUNT - 1;
			opColMax = 0;
		elseif site == 1
			mpodex = 1;
			opRowMax = 0;
			opColMax = OPCOUNT - 1;
		else
			mpodex = 2;
			opRowMax = OPCOUNT - 1;
			opColMax = OPCOUNT - 1;
		end

		B = mps{site};
		
		rowMax = size(B, 1);
		colMax = size(B, 2);

		conjB = zeros(colMax, rowMax, HILBY);
		for localState = 1 : 1 : HILBY
			conjB(:, :, localState) = ctranspose( B(:,:,localState) );
		end
	
		%rightBlock = sym( zeros(rowMax, OPCOUNT, rowMax) );		% SYM FOR DEBUG PURPOSES ONLY
		rightBlock = zeros(rowMax, OPCOUNT, rowMax);  

		for row = 1 : 1 : rowMax
			for opRow = 0 : 1 : opRowMax
				for conjCol = 1 : 1 : rowMax
					BWFB = 0;
					for ketState = 1 : 1 : HILBY
						for col = 1 : 1 : colMax
							WFB = 0;
							for braState = 1 : 1 : HILBY
								for opCol = 0 : 1 : opColMax
									FB = 0;
									for conjRow = 1 : 1 : colMax
										FB = FB + inner(conjRow, opCol + 1, col) * conjB(conjRow, conjCol, braState);
									end
									WFB = WFB + mpo{mpodex}(opRow * HILBY + ketState, opCol * HILBY + braState) * FB;
								end	% opCol
							end	% braState
							BWFB = BWFB + B(row, col, ketState) * WFB;
						end	% col
					end	% ketState
					rightBlock(row, opRow + 1, conjCol) = rightBlock(row, opRow + 1, conjCol) + BWFB;
				end 	% conjCol
			end	%opRow
		end	% row
		inner = rightBlock;
		fprintf('site %d contracted\n', site);
	end	% site  
end	% function
