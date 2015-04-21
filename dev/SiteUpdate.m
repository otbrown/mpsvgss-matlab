% SiteUpdate.m
% function which maps the eigenvector of the effective Hamiltonian on the target site back into a MPS matrix
% Oliver Thomson Brown
% 2015-03-31
%
% [RETURN]
% updateMatrix	: rowMax x colMax x HILBY double array, the rank-3 MPS tensor which minimises the energy on the target site
%
% [INPUTS]
% eigVec	: double array, the lowest eigenvalued eigenvector of the effective Hamiltonian on the target site, which is to be reshaped to be the new MPS tensor on that site
% rowMax	: double (integer), the maximum row index of the target MPS matrix
% colMax	: double (integer), the maximum column index of the target MPS matrix
% HILBY		: double (integer), the size of the local Hilbert space

function [ updateMatrix ] = SiteUpdate(eigVec, rowMax, colMax, HILBY)
	vecReshaped = reshape(eigVec, [colMax, rowMax, HILBY]);
	updateMatrix = permute(vecReshaped, [2, 1, 3]);
end
