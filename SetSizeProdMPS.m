% ProdMPS.m
% function to quickly construct an MPS representation of a product state
% Oliver Thomson Brown
% 20/10/2014
%
% [RETURN]
% matrices		: cell array, A_j_n = matrices{n}(:,:,j), where j indexes state, n indexes site
%
% [INPUT]
% stateArray		: column vector containing states i.e. [1; 0; 0; 1; 0; ...]
% HILBY			: dimension of Hilbert space on each site i.e. HILBY = 2 for a spin-half system

function [matrices] = ProdMPS(stateArray, HILBY, SIZE)

	% INPUT CHECK
	if size(stateArray,2) > 1
		disp('This function only accepts a column vector containing the state at each site in the chain.');
		return;
	end
	if SIZE < 2
		disp('Please enter a SIZE greater than or equal to 2');
		return;
	end

	stateArray = stateArray + 1;		% corrects indexing (states begin at 0, MATLAB begins at 1...)

	% CONSTANT GATHERING
	L = size(stateArray,1);
	
	% RETURN ALLOCATION
	matrices = cell(L,1);
	initMatrix = zeros(SIZE,SIZE,HILBY);
	[matrices{:}] = deal(initMatrix);
 
	aOne = eye(1,2);
	aL = eye(2,1);
	
	% CONSTRUCTION
	matrices{1}(1,1:2,stateArray(1)) = aOne;
	matrices{L}(1:2,1,stateArray(L)) = aL;

	rowSize = 2;
	if L == 3
		colSize = 2;
	else
		colSize = min(SIZE,4);
	end

	for site = 2 : 1 : L-1			% bulk sites treated here		
		lLen = site + 1;
		rLen = L - site - 1;
		len = min(lLen,rLen);
		state = stateArray(site);
	
		matrices{site}(1:rowSize,1:colSize,state) = eye(rowSize, colSize);		
		
		rowSize = colSize;
		colSize = min(SIZE, 2^len);
	end 	 

end
