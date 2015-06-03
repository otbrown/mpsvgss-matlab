% GrowLeft.m
% function which updates the contraction from the first site of a tensor network to include the next site to the right -- updates the leftBlock to include the next site along
% Oliver Thomson Brown
% 2015-04-10
%
% [Return]
% updateBlock	: 3-D double array, contains the new tensor contraction from the left
% 
% [Inputs]
% siteTensor	: 3-D double array, the MPS tensor on the site which is to be included in the left block
% mpoTensor	: 2-D double array, the MPO tensor for the site which is to be included in the left block
% leftBlock	: 3-D double array, the block which is to be 'grown' by including the contraction through the target site
% rowMax	: integer, the number of rows in the mps tensor on the site which is to be included in the left block
% colMax	: integer, the number of columns in the mps tensor on the site which is to be included in the right block
% HILBY		: integer, the size of the local state space
% opRowMax	: integer, the number of rows in the block notation form of the matrix product operator counting from 0
% opColMax	: integer, the number of columns in the block notation form of the matrix product opereator counting from 0
% OPCOUNT	: integer, the larger of opRowMax and opColMax plus one

function [ updateBlock ] = GrowLeft(siteTensor, mpoTensor, leftBlock, rowMax, colMax, HILBY, opRowMax, opColMax, OPCOUNT)
	% pre-allocate return array
	updateBlock = complex(zeros(colMax, colMax, OPCOUNT));
	       
	% perform contraction
	for opCol = 0 : 1 : opColMax
		for col = 1 : 1 : colMax
			AWFA = 0;
			for braState = 1 : 1 : HILBY
				for conjCol = 1 : 1 : rowMax
					WFA = 0;
					for ketState = 1 : 1 : HILBY
						for opRow = 0 : 1 : opRowMax
							WFA = WFA + mpoTensor(opRow * HILBY + braState, opCol * HILBY + ketState) * leftBlock(conjCol, :, opRow + 1) * siteTensor(:, col, ketState);
						end % opRow
					end % ketState
					AWFA = AWFA + conj( siteTensor(conjCol, 1 : colMax, braState) ) * WFA;
				end % conjCol
			end % braState
			updateBlock(1 : colMax, col, opCol + 1) = updateBlock(1 : colMax, col, opCol + 1) + transpose(AWFA);
		end % col
	end % opCol
end
