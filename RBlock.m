% RBlock.m
% Oliver Thomson Brown
% 2015-01-15
% DOCSTRING!

function [ rightBlock ] = RBlock(mps, mpo, TARGET)

	% PULL DATA LIKE TASHA YAR
	L = size(mps, 1);
	HILBY = size(mps{L}, 3);
	OPCOUNT = size(mpo{3}, 1) / HILBY;

	inner = ones(2, OPCOUNT, 2);

	for site = L : -1 : TARGET + 1
		disp(site);	% DEBUG
				
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
		conjB = cat(3, ctranspose( mps{site}(:, :, 1) ), ctranspose( mps{site}(:, :, 2) ) );

		rowMax = size(B, 1);
		colMax = size(B, 2);
	
		rightBlock = sym( zeros(rowMax, OPCOUNT, rowMax) );		% SYM FOR DEBUG PURPOSES ONLY
  
		for row = 1 : 1 : rowMax
			for opCol = 0 : 1 : opColMax
				for conjCol = 1 : 1 : rowMax
					BWFB = 0;
					for ketState = 1 : 1 : HILBY
						for col = 1 : 1 : colMax
							WFB = 0;
							for braState = 1 : 1 : HILBY
								for opRow = 0 : 1 : opRowMax
									FB = 0;
									for conjRow = 1 : 1 : colMax
										FB = FB + inner(conjRow, opRow + 1, col) * conjB(conjRow, conjCol, braState);
									end
									WFB = WFB + mpo{mpodex}(opRow * HILBY + ketState, opCol * HILBY + braState) * FB;
								end	% opRow
							end	% braState
							BWFB = BWFB + B(row, col, ketState) * WFB;
						end	% col
					end	% ketState
					rightBlock(row, opCol + 1, conjCol) = rightBlock(row, opCol + 1, conjCol) + BWFB;
				end 	% conjCol
			end	%opCol
		end	% row
		
		inner = rightBlock;
	end	% site  
end	% function
