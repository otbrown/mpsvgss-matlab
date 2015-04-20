% Decomp.m
% Function to filter decomposition of a state by SVD. If exact state is required full SVD is faster. 
% If SVD compression is required, sparse SVD is faster and allows compression to be applied. 
% Oliver Thomson Brown
% 16/10/2014

function [U, S, V] = Decomp(matrix, COMPRESS)

if COMPRESS == Inf				% COMPRESS == Inf corresponds to exact state decomposition
	[U, S, V] = svd(full(matrix));
else
	[U, S, V] = svds(matrix, COMPRESS);
end
