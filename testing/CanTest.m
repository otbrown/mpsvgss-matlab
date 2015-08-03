% CanTest.m
% a class to test the function Can.m or more specifically,
% its sub-functions LCan.m and RCan.m
% Oliver Thomson Brown
% 15-07-20

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) CanTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        testMPS;
        testMPSLength;
        testHILBY;
    end

    properties (MethodSetupParameter)
        HILBY = {2, 2, 3, 4, 5, 2, 3, 4, 5};
        mpsLength = {5, 7, 6, 5, 5, 7, 6, 5, 5};
        COMPRESS = {0, 0, 0, 0, 0, 2, 6, 12, 15};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(testCase, HILBY, mpsLength, COMPRESS)
            % create mps
            testCase.testMPS = CompMPS(HILBY, mpsLength, COMPRESS);
            testCase.testMPSLength = mpsLength;
            testCase.testHILBY = HILBY;
        end
    end

    methods (Test)
        %function testLNorm(testCase)
        %    % bring in to left-canonical form
        %    testCase.testMPS = Can(testCase.testMPS, 1 : 1 : (testCase.testMPSLength - 1), 'L');
        %    % for each site and local state, assert A'*A = I
        %    for site = 1 : 1 : (testCase.testMPSLength - 1)
        %        matrix = testCase.testMPS{site};
        %        for localState = 1 : 1 : testCase.testHILBY
        %            colSz = size(matrix, 2);
        %            hermProd = ctranspose(matrix(:, :, localState)) * matrix(:, :, localState);
        %            testCase.assertEqual(hermProd, eye(colSz), 'AbsTol', testCase.absTol);
        %        end
        %    end
        %end

        %function testRNorm(testCase)
        %    % bring in to right-canonical form
        %    testCase.testMPS = Can(testCase.testMPS, testCase.testMPSLength : -1 : 2, 'R');
        %    % for each site and local state assert A*A' = I
        %    for site = testCase.testMPSLength : -1 : 2
        %        matrix = testCase.testMPS{site};
        %        for localState = 1 : 1 : testCase.testHILBY
        %            rowSz = size(matrix, 1);
        %            hermProd = matrix(:, :, localState) * ctranspose(matrix(:, :, localState));
        %            testCase.assertEqual(hermProd, eye(rowSz), 'AbsTol', testCase.absTol);
        %        end
        %    end
        %end

        function testNormL(testCase)
            % bring in to left-canonical form
            testCase.testMPS = Can(testCase.testMPS, 1 : 1 : (testCase.testMPSLength - 1), 'L');
            % rebuild state-vector and assert psi' * psi = 1
            stateVec = Rebuild(testCase.testMPS);
            norm = ctranspose(stateVec) * stateVec;
            testCase.assertEqual(norm, 1, 'AbsTol', testCase.absTol);
        end

        function testNormR(testCase)
            % bring in to right-canonical form
            testCase.testMPS = Can(testCase.testMPS, testCase.testMPSLength : -1 : 2, 'R');
            % rebuild state-vector and assert norm equal to 1
            stateVec = Rebuild(testCase.testMPS);
            norm = ctranspose(stateVec) * stateVec;
            testCase.assertEqual(norm, 1, 'AbsTol', testCase.absTol);
        end
    end

end
