% RCan.m
% function which right-canonises an already formed matrix product state
% please note that the function Can acts as an interface to this one
% in particular Can ensures that the supplied route is valid 
% Oliver Thomson Brown
% 13/11/2014
%
% [RETURN]
% rmps		: cell array containing right-canonised MPS, A_state_site = lmps{site}(:,:,state)
%
% [INPUTS]
% mps		: cell array containing a matrix product state, A_state_site = mps{site}(:,:,state)
% route		: vector containing sites which are to be canonised, if whole chain is to be canonised 
%		  route should be L : -1 : 2 where L is the length of the chain

function [rmps] = RCan(mps, route)
	% RETURN ALLOCATION
	rmps = mps;
	
	for site = route		
		M = cat(2, rmps{site}(:,:,1), rmps{site}(:,:,2) );
		[U, S, V] = svd(M);

		V = ctranspose(V);

		Vcol = size(V,2);
		rowLim = size(M,1);
		colLim = size(M,2) / 2;

		rmps{site}(:, :, :) = cat(3, V(1:rowLim, 1:colLim), V(1:rowLim, Vcol/2 + 1 : Vcol/2 + colLim));

		chain = U * S;

		colLim = rowLim;
		rowLim = size(rmps{site-1}(:,:,1), 1);

		N = cat(3, rmps{site-1}(:,:,1) * chain, rmps{site-1}(:,:,2) * chain);
		
		rmps{site-1} = cat(3, N(1:rowLim,1:colLim,1), N(1:rowLim,1:colLim,2));

		%fprintf('site %d right-normalised\n', site);
	end

end
