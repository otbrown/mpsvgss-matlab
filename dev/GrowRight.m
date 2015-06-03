% GrowRight.m
% function which grows a contraction through from the right side of the system by including the contraction through a target site
% Oliver Thomson Brown
% 2015-04-14
%
% [RETURN]
% updateBlock	: 3-D double array, the contraction from the right hand side of the tensor network up to and including the TARGET site
%
% [INPUTS]
% siteTensor	: 3-D double array, the mps tensor on the TARGET site
% mpoTensor	: 2-D double array, the matrix product operator for the TARGET site
% rightBlock	: 3-D double array, the contraction from the last site through the tensor network up to but *not* including the TARGET site
% rowMax	: integer, the number of rows in siteTensor
% colMax	: integer, the number of columns in siteTensor
% HILBY		: integer, the dimension of the local state space
% opRowMax	: integer, the number of rows in the block representation of the matrix product operator, counting from 0
% opColMax	: integer, the number of columns in the block representation of the matrix product operator, counting from 0
% OPCOUNT	: integer, the larger of opRowMax and opColMax plus one

function [ updateBlock ] = GrowRight(siteTensor, mpoTensor, rightBlock, rowMax, colMax, HILBY, opRowMax, opColMax, OPCOUNT)
	% pre-allocate return array
	updateBlock = complex(zeros(rowMax, rowMax, OPCOUNT));
	
	% transpose siteTensor to improve inner-loop efficiency
    	tSiteTensor = permute(siteTensor, [2, 1, 3]);

	% perform contraction
	for opRow = 0 : 1 : opRowMax
		for row = 1 : 1 : rowMax
			BWFB = 0;
			for braState = 1 : 1 : HILBY
				for conjRow = 1 : 1 : colMax
					WFB = 0;
					for ketState = 1 : 1 : HILBY
						for opCol = 0 : 1 : opColMax
							WFB = WFB + mpoTensor(opRow * HILBY + braState, opCol * HILBY + ketState) * rightBlock(conjRow, :, opCol + 1) * tSiteTensor(:, row, ketState);
						end % opCol
					end % ketState 
					BWFB = BWFB + conj(tSiteTensor(conjRow, 1 : rowMax, braState)) * WFB;
				end % conjRow
			end % braState
			updateBlock(1 : rowMax, row, opRow + 1) = updateBlock(1 : rowMax, row, opRow + 1) + transpose(BWFB);
		end % row
	end % opRow
end
