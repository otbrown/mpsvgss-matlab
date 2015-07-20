% CanTest.m
% a class to test the function Can.m or more specifically,
% its sub-functions LCan.m and RCan.m
% Oliver Thomson Brown
% 15-07-20

classdef CanTest < matlab.unittest.TestCase

    properties (TestParameter)
        HILBY = [2, 2, 3, 4, 5, 2, 3, 4, 5];
        mpsLength = [5, 7, 6, 5, 5, 7, 6, 5, 5];
        COMPRESS = [0, 0, 0, 0, 0, 2, 6, 12, 15]; 
    end

    methods (Test, ParameterCombination='sequential')
        function testLNorm(testCase, HILBY, mpsLength, COMPRESS)
            absTol = 1E-14;
            % create mps
            mps = CompMPS(HILBY, mpsLength, COMPRESS);
            % bring in to left-canonical form
            mps = Can(mps, 1 : 1 : (mpsLength - 1), 'L');
            % for each site and local state, assert A'*A = I
            for site = 1 : 1 : (mpsLength - 1)
                matrix = mps{site};
                for localState = 1 : 1 : HILBY
                    colSz = size(matrix, 2);
                    hermProd = ctranspose(matrix) * matrix;
                    testCase.assertEqual(hermProd, eye(colSz), 'AbsTol', absTol);
                end
            end
        end
    end

end   
