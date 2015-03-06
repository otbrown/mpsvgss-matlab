% LCan.m
% function which left-canonises an already formed matrix product state
% please note that the function Can acts as an interface to this one
% in particular Can ensures that the supplied route is valid 
% Oliver Thomson Brown
% 10/11/2014
%
% [RETURN]
% lmps		: cell array containing left-canonised MPS, A_state_site = lmps{site}(:,:,state)
%
% [INPUTS]
% mps		: cell array containing a matrix product state, A_state_site = mps{site}(:,:,state)
% route		: vector containing sites which are to be canonised, if whole chain is to be canonised 
%		  route should be 1 : 1 : L-1 where L is the length of the chain

function [lmps] = LCan(mps, route)

	threshold = 1E-14;
	
	% RETURN ALLOCATION
	lmps = mps;
	
	for site = route		
		M = cat(1, lmps{site}(:,:,1), lmps{site}(:,:,2) );
		[U, S, V] = svd(M);

		U(abs(U) < threshold) = 0;
		S(abs(S) < threshold) = 0;
		V(abs(V) < threshold) = 0;

		Urow = size(U,1);
		rowLim = size(M,1) / 2;
		colLim = size(M,2);

		lmps{site}(:, :, :) = cat(3, U(1:rowLim, 1:colLim), U(Urow/2 + 1 : Urow/2 + rowLim, 1:colLim));

		chain = S * ctranspose(V);

		rowLim = colLim;
		colLim = size(lmps{site+1}(:,:,1), 2);

		N = cat(3, chain * lmps{site+1}(:,:,1), chain * lmps{site+1}(:,:,2));
		
		lmps{site+1} = cat(3, N(1:rowLim,1:colLim,1), N(1:rowLim,1:colLim,2));

		%fprintf('site %d left-normalised\n', site);
	end

end
