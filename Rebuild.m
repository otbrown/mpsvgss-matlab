% Rebuild.m
% function which rebuilds the explicit (sparse) state vector from the MPS form -- HILBY = 2 ONLY
% for debugging use only -- this will never be efficient! 
% Oliver Thomson Brown
% 16/10/2014
%
% [RETURNS]
% state		: d^L double array, where d is the hilbert space of a single site and L is the number of sites, COMPULSORY
% epsilon	: d^L double array, contains the difference between the 'real' state vector and the rebuilt one OPTIONAL
%
% [INPUTS]
% matrices	: L by 1 cell array, contains the MPS coefficient matrices
% init		: d^L double array, contains the initial state vector which was decomposed

function [state, varargout] = Rebuild(matrices, varargin)

	% INPUT/RETURN CHECKS & TIDY-UP 
	nargoutchk(0,2);	
	narginchk(1,2);			% checks number of arguments for input and output

	threshold = 1E-14;
	L = size(matrices, 1);
	SPACE = 2^L;
	
	if nargout == 2
		if nargin == 1
			disp('You must provide the original state vector if you wish to create the error vector, epsilon.');
			return;
		else
			init = varargin{1};	% varargin{1} is the user-provided original state vector -- normalisation check follows
			normCheck = ctranspose(init) * init;
			if normCheck ~= 1
				init = init / sqrt(normCheck);
			end
		end
		epsilon = sparse(SPACE,1);
	end

	state = sparse(SPACE, 1);
			
	% FORM BASIS STATES AND BINARY COMBINATIONS -- SCALES HORRIBLY
	bins = zeros( SPACE, L, 'uint8' );	% SPACE by L uint8 array, will contain binary combos
	comp = speye( SPACE );       		% computational basis states in sparse form

	for sigma = 1 : SPACE			% form all possible binary combinations
	    bins( sigma, : ) = double( dec2bin( sigma-1, L ) - 48 );
	end
	bins = flipud(bins);    		% binary ordering matches colwise comp states

	% REBUILD STATE VECTOR -- SCALES HORRIBLYER	
	for sigma = 1 : SPACE
	    coefft = 1 ;
	    
	    for site = 1 : 1 : L
		coefft = coefft*matrices{site}(:, :, 2^bins(sigma , site));
	    end
	    
	    state = state + coefft*comp(:, sigma);
	end

	state(abs(state) < threshold) = 0;

	if nargout == 2
		epsilon = abs(state - init);
		epsilon(epsilon < threshold) = 0;
		varargout{1} = epsilon;		% varargout{1} is the difference between the rebuilt and original states 
	end

end
