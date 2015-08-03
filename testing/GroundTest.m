% GroundTest.m
% testing class to check top-level ground-state search function
% IMPORTANT: REQUIRES FILE CONTAINING THE MPO YOU WISH TO TEST AGAINST
% --also the given stateVec should match that MPO...
% Oliver Thomson Brown
% 2015-08-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) GroundTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        stateVecTol = 1E-8;
        testMPO = load('MPO_on-site+hop.mat', 'mpo');
        testMPS;
        testGroundMPS;
        testHILBY;
        testMPSLength;
    end

    properties (MethodSetupParameter)
        HILBY = {2, 2, 2, 2, 2, 2}
        mpsLength = {5, 7, 9, 7, 9, 12}
        COMPRESS = {0, 0, 0, 4, 8, 16}
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(testCase, HILBY, mpsLength, COMPRESS)
            % create mps
            testCase.testMPS = CompMPS(HILBY, mpsLength, COMPRESS);
            [testCase.testGroundMPS, ~] = Ground(testCase.testMPS, testCase.testMPO.mpo, 1E-12, 100);
            testCase.testHILBY = HILBY;
            testCase.testMPSLength = mpsLength;
        end
    end

    methods (TestMethodTeardown)
        function cleanCheckpointFiles(testCase)
            delete('*chkpnt.mat');
        end
    end

    methods (Test)
        function testClass(testCase)
            testCase.assertClass(testCase.testGroundMPS, 'cell');
        end

        function testCellShape(testCase)
            testCase.assertSize(testCase.testGroundMPS, [testCase.testMPSLength, 1]);
        end

        function testLocalStateSpace(testCase)
            testCase.assertEqual(size(testCase.testGroundMPS{1}, 3), testCase.testHILBY);
        end

        function testGroundNorm(testCase)
            stateVec = Rebuild(testCase.testGroundMPS);
            norm = ctranspose(stateVec) * stateVec;
            testCase.assertEqual(norm, 1, 'AbsTol', testCase.absTol);
        end

        function testGroundState(testCase)
            stateVec = Rebuild(testCase.testGroundMPS);
            SPACE = testCase.testHILBY ^ testCase.testMPSLength;
            testCase.assertEqual(abs(stateVec), eye(SPACE, 1), 'AbsTol', testCase.stateVecTol);
        end
    end
end
