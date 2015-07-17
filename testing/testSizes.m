% testSizes.m
% function to check that the matrices in a supplied matrix product state are sized consistently
% Oliver Thomson Brown
% 2015-07-17
%
% Inputs:
% testMPS   : cell array, the mps to be tested
% HILBY     : the correct local Hilbert space dimension
% COMPRESS  : the maximum dimension of the matrices, if 0 is supplied it is converted to infinity -- an exact MPS

function testSizes(testMPS, HILBY, COMPRESS)
    if COMPRESS == 0
        COMPRESS = Inf;
    end
    % initialise rowSz and colSz
    rowSz = 1;
    colSz = min(HILBY, COMPRESS);
    % sweep through each site in the matrix product state
    for site = 1 : 1 : length(testMPS) - 1;    
        matSz = size(testMPS{site});
        assert(matSz(1) == rowSz, 'Error: Matrices are inconsistently sized.');
        assert(matSz(2) == colSz, 'Error: Matrices are inconsistently sized.');
        assert(matSz(3) == HILBY, 'Error: Matrices are inconsistently sized.');
        % gauge whether there are more sites to the left, or to the right of the current one
        lLen = site + 1;
        rLen = length(testMPS) - site - 1;
        len = min(lLen, rLen);
        % set next site rowSz and colSz
        rowSz = colSz;
        colSz = min(COMPRESS, HILBY^len);
    end
    matSz = size(testMPS{end});
    assert(matSz(1) == min(HILBY, COMPRESS), 'Error: Matrices are inconsistently sized.');
    assert(matSz(2) == 1, 'Error: Matrices are inconsistently sized.');
    assert(matSz(3) == HILBY, 'Error: Matrices are inconsistently sized.');
end
