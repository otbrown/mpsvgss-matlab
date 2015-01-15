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
	
		rightBlock = sym( zeros(colMax, OPCOUNT, colMax) );		% SYM FOR DEBUG PURPOSES ONLY
  
		for conjRow = 1 : 1 : colMax
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
								end	% opRow
							end	% ketState
							BWFB = BWFB + conjB(conjRow, opCol + 1, col) * WFB;
						end	% conjCol
					end	% braState
					rightBlock(conjRow, opCol + 1, col) = rightBlock(conjRow, opCol + 1, col) + BWFB;
				end 	% col
			end	%opCol
		end	%conjRow
		
		inner = rightBlock
	end	% site  
end	% function
