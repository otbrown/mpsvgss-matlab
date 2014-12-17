% Expect.m
% a function to calculate the expectation value of a matrix product state
% contracts alternately through physical indices chain then virtual indices at each site
% Oliver Thomson Brown
% 21/10/2014

% [RETURN]
% expectation		: the result of the calculation, this may either be a scalar, vector, or matrix, depending on the route through the chain

% [INPUT]
% braMPS		: cell array containing the conjugate MPS coefficient matrices
% ketMPS		: cell array containing the MPS coefficient matrices
% operator		: cell array containing the operator at each site
% route			: row vector containing the route through the chain, i.e. route = 1:1:L

function expectation = Expect(braMPS, ketMPS, operator, route)
	threshold = 1E-14; 		% numerical error threshold

	% CONSTANT GATHERING
	L = size(ketMPS,1);
	HILBY = size(ketMPS{1},3);

	% CONTRACTION
	inner = 1;			% <braMPS|ketMPS> =  M'_L ... M'_1 * M_1 ... M_L so code contracts in this manner, keeping track of what's in the middle!
	for site = route;
		outer = 0;
		for conjState = 1:1:HILBY
			for state = 1:1:HILBY
				opRow = HILBY + 1 - conjState;
				opCol = HILBY + 1 - state;
				outer = outer + operator{site}(opRow,opCol) * ctranspose( braMPS{site}(:,:,conjState) ) * inner * ketMPS{site}(:,:,state);
			end
		end
		outer(abs(outer) < threshold) = 0;		% stability measure
		inner = outer;
	end
	expectation = outer;
end
