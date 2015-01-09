% LBlock.m
% function which
% Oliver Thomson Brown
% 2015-01-08
% 
%	[RETURN]
%	leftBlock	:	double array, contains contraction from first site through to site left of target
%
%	[INPUTS]
%	mps			:	cell array, contains matrix product state
%	mpo			:	cell array, contains matrix product operator mpo{1} is first site op, mpo{2} is bulk site op, mpo{3} is last site op
%	TARGET		:	constant int, contains target site	

function [ leftBlock ] = LBlock(mps, mpo, TARGET)

	% pull data from inputs
	L = size(mps, 1);
	HILBY = size(mps{1}, 3);
	OPCOUNT = size(mpo{1}, 2) / HILBY;
	
	inner = ones(1, OPCOUNT, 1);	% initialise '0th site'

	for site = 1 : 1 : TARGET - 1
		disp(site);		%DEBUG!
		if site == 1		% select correct MPO matrix
			mpodex = 1;
			opRowMax = 0;
			opColMax = OPCOUNT - 1;
		elseif site == L
			mpodex = 3;
			opRowMax = OPCOUNT -1;
			opColMax = 0;
		else
			mpodex = 2;
			opRowMax = OPCOUNT - 1;
			opColMax = OPCOUNT - 1;
		end

		A = mps{site};
		conjA = cat(3, ctranspose(mps{site}(:,:,1)), ctranspose(mps{site}(:,:,2)));
		
		rowMax = size(A, 1);
		colMax = size(A, 2);

		leftBlock = sym(zeros(colMax, OPCOUNT, colMax));	% REMOVE SYM WHEN FINISHED -- DEBUG!

		for conjRow = 1 : 1 : colMax				% SUMMATIONS ORDERED AS PRESENTED IN SCHOLLWOECK p65
			for opCol = 0 : 1 : opColMax
				for col = 1 : 1 : colMax
					AWFA = 0;
					for braState = 1 : 1 : HILBY
						for conjCol = 1 : 1 : rowMax
							WFA = 0;
							for ketState = 1 : 1 : HILBY
								for opRow = 0 : 1 : opRowMax
									FA = 0;
									for row = 1 : 1 : rowMax
										FA = FA + inner(conjCol, opRow + 1, row) * A(row, col, ketState);
									end
									WFA = WFA + mpo{mpodex}(opRow * HILBY + braState, opCol * HILBY + ketState) * FA;
								end
							end
							AWFA = AWFA + conjA(conjRow, conjCol, braState) * WFA;
						end		% conjCol
					end 		% braState
					leftBlock(conjRow, opCol + 1, col) = leftBlock(conjRow, opCol + 1, col) + AWFA;
				end		% col
			end		%opCol
		end		% conjRow
		inner = leftBlock;
	end		% site

end		% function
