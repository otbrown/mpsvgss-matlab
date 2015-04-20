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
	updateMatrix = zeros(rowMax, colMax, HILBY);

	for ketState = 0 : 1 : HILBY - 1
		for row = 0 : 1 : rowMax - 1
			for col = 1 : 1 : colMax
				vecdex = ketState * rowMax * colMax + row * colMax + col;
				updateMatrix(row + 1, col, ketState + 1) = eigVec(vecdex);
			end
		end
	end
end
