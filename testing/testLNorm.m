% testLNorm.m
% function to assert the supplied testMPS is left-normalised
% Oliver Thomson Brown
% 2015-07-20
%
% Inputs:
% testCase  : an interactive test-case, object tc = matlab.unittest.TestCase.forInteractiveUse;
% testMPS   : cell array containing the matrix product state to be tested
% tolerance : accuracy to which arrays should be compared

function testLNorm(testCase, testMPS, tolerance);
    mpsLength = size(testMPS, 1);
    for site = 1 : 1 : mpsLength - 1
        [~, colSz, HILBY] = size(testMPS{site});
        for localState = 1 : 1 : HILBY
            matrix = testMPS{site}(:, :, localState);
            hermProd = ctranspose(matrix) * matrix;
            assertEqual(testCase, hermProd, eye(colSz), 'AbsTol', tolerance);
        end
    end
end
