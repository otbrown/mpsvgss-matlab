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

	if TARGET == 0
		inner = ones(1, OPCOUNT, 1);
	else
		targRow = size( mps{TARGET}, 1 );
		targCol = size( mps{TARGET}, 2 );
		inner = ones(targRow, OPCOUNT, targRow);
	end

	for site = TARGET + 1 : 1 : L
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

		fprintf('site %d rowMax = %d\n', site, rowMax);
		fprintf('site %d colMax = %d\n', site, colMax);

		conjB = zeros(colMax, rowMax, HILBY);
		for localState = 1 : 1 : HILBY
			conjB(:, :, localState) = ctranspose( B(:,:,localState) );
		end
	
		%rightBlock = sym( zeros(rowMax, OPCOUNT, rowMax) );		% SYM FOR DEBUG PURPOSES ONLY
		rightBlock = zeros(colMax, OPCOUNT, colMax);  
%{
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
%}

		for conjRow = 1 : 1 : colMax				% SUMMATIONS ORDERED AS PRESENTED IN SCHOLLWOECK p65
			for opCol = 0 : 1 : opColMax
				for col = 1 : 1 : colMax
					BWFB = 0;
					for braState = 1 : 1 : HILBY
						for conjCol = 1 : 1 : rowMax
							WFB = 0;
							for ketState = 1 : 1 : HILBY
								for opRow = 0 : 1 : opRowMax
									FB = 0;
									for row = 1 : 1 : rowMax
										FB = FB + inner(conjCol, opRow + 1, row) * B(row, col, ketState);
									end
									WFB = WFB + mpo{mpodex}(opRow * HILBY + braState, opCol * HILBY + ketState) * FB;
								end
							end
							BWFB = BWFB + conjB(conjRow, conjCol, braState) * WFB;
						end		% conjCol
					end 		% braState
					rightBlock(conjRow, opCol + 1, col) = rightBlock(conjRow, opCol + 1, col) + BWFB;
				end		% col
			end		%opCol
		end		% conjRow
		inner = rightBlock;
		fprintf('site %d contracted\n', site);
	end	% site  
end	% function
