% RBlock.m
% Oliver Thomson Brown
% 2015-02-10
%
% [RETURN]
% rightBlock	: rank 3 tensor, contains contraction of tensor network from last site up to the site to the TARGET site
%
% [INPUTS]
% mps		: cell array, matrix product state, mps{site} contains the rank 3 site tensor 
% mpo		: cell array, matrix product operator, mpo{1} coontains the first-site mpo, mpo{2} contains the bulk site mpo, mpo{3} contains the last site mpo
% TARGET	: int, location of the site up to which the right-block should be formed -- typically the next site to be updated

function [ rightBlock ] = RBlock( mps, mpo, TARGET )
	% pull data from inputs
	L = size( mps, 1 );

	if TARGET == L		% protection
		rightBlock = 1;
	else
		HILBY = size( mps{L}, 3);
		OPCOUNT = size( mpo{3}, 1) / HILBY;

		for site = L : -1 : TARGET + 1
			if site == L
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

			rowMax = size( B, 1);
			colMax = size( B, 2);
			
			conjB = zeros( colMax, rowMax, HILBY);
			for localState = 1 : 1 : HILBY
				conjB( :, :, localState) = ctranspose( B( :, :, localState) );
			end

			rightBlock = zeros( rowMax, rowMax, OPCOUNT);

			for conjCol = 1 : 1 : rowMax
				for opRow = 0 : 1 : opRowMax
					for row = 1 : 1 : rowMax
						BWFB = 0;
						for braState = 1 : 1 : HILBY
							for conjRow = 1 : 1 : colMax
								WFB = 0;
								for ketState = 1 : 1 : HILBY
									for opCol = 0 : 1 : opColMax
										FB = 0;
										for col = 1 : 1 : colMax
											if site == L
												FB = FB + B(row, col, ketState);
											else
												FB = FB + inner(conjRow, col, opCol + 1) * B(row, col, ketState);
											end
										end
										WFB = WFB + mpo{mpodex}(opRow * HILBY + braState, opCol * HILBY + ketState) * FB;
									end
								end
								BWFB = BWFB + conjB(conjRow, conjCol, braState) * WFB;
							end
						end
						rightBlock(conjCol, row, opRow + 1) = rightBlock(conjCol, row, opRow + 1) + BWFB;
					end
				end
			end
			inner = rightBlock;
			%fprintf('site %d contracted\n', site);
		end	%site
	end	% if in case TARGET == L
end	%function
