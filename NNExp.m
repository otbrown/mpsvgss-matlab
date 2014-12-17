% NNExp.m
% function to return the expectation value of a nearest neighbour Operator
% (such as a Hamiltonian) on close-bounded one-dimensional lattice 
% Oliver Thomson Brown
% 14/11/2014
%
% [RETURN]
% expectation		: double, the expectation value of the operator
%
% [INPUT]
% HILBY			: double, the dimension of single-site Hilbert space
% siteOp		: the single-site operator, should be a HILBY x HILBY array
% interOp		: the interaction poerator, should be HILBY x HILBY x 2 array, contains the site i and site j operator concatenated along dim3
% braMPS		: the conjugate (bra) matrix product state		< bra | OPERATOR | ket >
% ketMPS		: the ket matrix product state

[ expectation ] = NNExp(HILBY, siteOp, interOp, braMPS, ketMPS)
	L = size(ketMPS,1);			% constant gathering
	chain = 1 : 1 : L;			% allocations
	route = chain;	
	expectation = 0;	
	I = eye(HILBY);				
	operatorString = cell(L, 1);
	[operatorString{:}] = deal(I);

	for site = route			% on-site terms
		operatorString{site} =  siteOp;								% set
		expectation = expectation + Expect(braMPS, ketMPS, operatorString, chain);		% calc
		operatorSring{site} = I;								% reset
	end

	route = 1 : 1 : L - 1;			% nearest neighbour interaction terms
	for site = route
		operatorString{site} = interOp(:,:,1);							% set
		operatorString{site + 1} = interOp(:,:,2);
		expectation = expectation + Expect(braMPS, ketMPS, operatorString, chain);		% calc 
		
		operatorString{site} = interOp(:,:,2);							% switch
		operatorString{site + 1} = interOp(:,:,1);
		expectation = expectation + Expect(braMPS, ketMPS, operatorString, chain);		% calc 
		
		operatorString{site} = I;								% reset
		operatorString{site + 1} = I;						
	end
end
