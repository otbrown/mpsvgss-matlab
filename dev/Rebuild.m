% Rebuild.m
% function to rebuild a state-vector from a matrix product state
% !IMPORTANT NOTE!: for debugging use only -- this will *never* be efficient! 
% Oliver Thomson Brown
% 16/10/2014
%
% [RETURNS]
% state		: d^L double array, where d is the hilbert space of a single site and L is the number of sites, COMPULSORY
% epsilon	: d^L double array, contains the difference between the 'real' state vector and the rebuilt one OPTIONAL
%
% [INPUTS]
% mps	: L by 1 cell array, contains the MPS coefficient matrices
% init		: d^L double array, contains the initial state vector which was decomposed

function [state, varargout] = Rebuild(mps, varargin)
	% INPUT/RETURN CHECKS & TIDY-UP 
	nargoutchk(0,2);	
	narginchk(1,2);			% checks number of arguments for input and output

	L = size(mps, 1);
    HILBY = size(mps{1}, 3);
	SPACE = HILBY^L;
	
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
		epsilon = complex(zeros(SPACE,1));
	end

	state = complex(zeros(SPACE, 1));

    % set coefft = 1, and create storage for selected matrices
    coefft = 1;
    stateMat = cell(L, 1);

    % create all 0's coefficient, since it's base HILBY representation is obvious
    for site = 1 : 1 :  L
        stateMat{site} = mps{site}(:, :, 1);
    end
    for site = 1 : 1 : L
        coefft = coefft * stateMat{site};
    end
    state(1) = coefft;
    
    % REBUILD STATE VECTOR
    for stateDex = 2 : 1 : SPACE
        coefft = 1;
        % generate string which is the base HILBY decomposition of stateDex
        localStr = dec2base(stateDex - 1, HILBY, L);
        % locate correct matrices from mps
        for site = 1 : 1 : L
            stateMat{site} = mps{site}(:, :, str2num(localStr(site)) + 1 );
        end
        % multiply them to obtain the coefficient
        for site = 1 : 1 : L
            coefft = coefft * stateMat{site};
        end
        state(stateDex) = coefft;
    end
    
	if nargout == 2
		epsilon = abs(state - init);
		varargout{1} = epsilon;		% varargout{1} is the difference between the rebuilt and original states 
	end        
end
