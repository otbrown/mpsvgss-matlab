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

function [ updateBlock ] = GrowLeft(siteTensor, mpoTensor, leftBlock)
	% gather constants
	[ rowMax, colMax, HILBY ] = size( siteTensor );
	opRowMax = ( size(mpoTensor, 1) / HILBY ) - 1;
	opColMax = ( size(mpoTensor, 2) / HILBY ) - 1;
	OPCOUNT = max(opRowMax, opColMax) + 1;
	leftBlockSize = size( leftBlock );

	% pre-allocate return array
	updateBlock = zeros(colMax, colMax, OPCOUNT);
	
	% perform nested sum
	for conjRow = 1 : 1 : colMax
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
									if leftBlockSize == 1
										FA = FA + siteTensor(row, col, ketState);
									else
										FA = FA + leftBlock(conjCol, row, opRow + 1) * siteTensor(row, col, ketState);
									end
								end % row
								WFA = WFA + mpoTensor(opRow * HILBY + braState, opCol * HILBY + ketState) * FA;
							end % opRow
						end % ketState
						AWFA = AWFA + conj( siteTensor(conjCol, conjRow, braState) ) * WFA;
					end % conjCol
				end % braState
				updateBlock(conjRow, col, opCol + 1) = updateBlock(conjRow, col, opCol + 1) + AWFA;
			end % col
		end % opCol
	end % conjRow
end
		 
