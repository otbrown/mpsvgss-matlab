% GrowBlock.m
% function which can be used to update the left or right contraction up to a given site by including a target site in the contraction up to it
% note that the TARGET supplied should be the site which is to be included, so for example the intended use of this function would be: 
% left{TARGET + 1} = GrowBlock(mps, mpo, TARGET, 'L', left)
% Oliver Thomson Brown
% 2015-04-14
%
% [RETURN]
% updateBlock	: 3-D double array, contraction through the whole tensor network from the first or last site up to and including the TARGET site
% 
% [INPUTS]
% mps		: cell array, contains the most recent matrix product state for the system
% mpo		: cell array, contains the matrix product operator tensors for the system
% TARGET	: integer, the index of the site which is to be included in the contraction
% SIDE		: 1-character string, 'L' or 'R', the side of the network which is being contracted through 'L' is from the first site, 'R' is from the last site
% blockCell	: cell array, each cell contains the contraction through the network up to that site so left{1} for example contains 1, since there are no sites left of the first

function [ updateBlock ] = GrowBlock(mps, mpo, TARGET, SIDE, blockCell)
	% gather mps data
	L = length(mps);
	siteTensor = mps{TARGET};
	[rowMax, colMax, HILBY] = size(siteTensor);
	
	% select correct mpo tensor
	if TARGET == 1
		mpodex = 1;
	elseif TARGET == L
		mpodex = 3;
	else
		mpodex = 2;
	end
	mpoTensor = mpo{mpodex};

	% gather mpo data
	opRowMax = ( size(mpoTensor, 1) / HILBY ) - 1;
	opColMax = ( size(mpoTensor, 2) / HILBY ) - 1;
	OPCOUNT = max(opRowMax, opColMax) + 1;

	% select block
	block = blockCell{TARGET};
	noBlock = 0;

	% make Grow call
	if SIDE == 'L'
		if TARGET == 1
			noBlock = 1;
		end
		updateBlock = GrowLeft(siteTensor, mpoTensor, block, noBlock, rowMax, colMax, HILBY, opRowMax, opColMax, OPCOUNT);
	elseif SIDE == 'R'
		if TARGET == L
			noBlock = 1;
		end
		updateBlock = GrowRight(siteTensor, mpoTensor, block, noBlock, rowMax, colMax, HILBY, opRowMax, opColMax, OPCOUNT);
	else
		updateBlock = 0;
		fprintf('Please enter a valid SIDE of the network (L or R).\n')
	end
end  
