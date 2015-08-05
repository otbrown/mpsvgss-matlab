% CanTest.m
% a class to test the function Can.m or more specifically,
% its sub-functions LCan.m and RCan.m
% Oliver Thomson Brown
% 2015-07-20

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

        function testNormL(testCase)
            % bring in to left-canonical form
            testCase.testMPS = Can(testCase.testMPS, 1 : 1 : (testCase.testMPSLength - 1), 'L');
            % build identity mpo
            impo = [{eye(testCase.testHILBY)};{eye(testCase.testHILBY)};{eye(testCase.testHILBY)}];
            % build left contraction to last site
            lBlock = cell(testCase.testMPSLength, 1);
            lBlock{1} = 1;
            for site = 2 : testCase.testMPSLength
                lBlock{site} = GrowBlock(testCase.testMPS, impo, ...
                    (site - 1), 'L', lBlock);
            end
            % calculate norm of mps
            norm = Expect(testCase.testMPS, impo, lBlock{testCase.testMPSLength}, 1, testCase.testMPSLength);
            % assert norm == 1
            testCase.assertEqual(norm, 1, 'AbsTol', testCase.absTol);
        end
        
        function testStateL(testCase)
            % rebuild state vector
            stateVec_init = Rebuild(testCase.testMPS);
            % bring in to left-canonical form
            testCase.testMPS = Can(testCase.testMPS, 1 : 1 : (testCase.testMPSLength - 1), 'L');
            % build new state vector
            stateVec_L = Rebuild(testCase.testMPS);
            % assert they are the same up to some phase (imaginary
            % components may differ)
            testCase.assertEqual(abs(stateVec_init), abs(stateVec_L), 'AbsTol', testCase.absTol);
        end

        function testNormR(testCase)
            % bring in to right-canonical form
            testCase.testMPS = Can(testCase.testMPS, testCase.testMPSLength : -1 : 2, 'R');
            % build identity mpo
            impo = [{eye(testCase.testHILBY)};{eye(testCase.testHILBY)};{eye(testCase.testHILBY)}];
            % build right contraction to first site
            rBlock = cell(testCase.testMPSLength, 1);
            rBlock{testCase.testMPSLength} = 1;
            for site = testCase.testMPSLength - 1 : -1 : 1
                rBlock{site} = GrowBlock(testCase.testMPS, impo, ...
                    (site + 1), 'R', rBlock);
            end
            % calculate norm of mps
            norm = Expect(testCase.testMPS, impo, 1, rBlock{1}, 1);
            % assert norm == 1
            testCase.assertEqual(norm, 1, 'AbsTol', testCase.absTol);
        end
        
        function testStateR(testCase)
            % rebuild state vector
            stateVec_init = Rebuild(testCase.testMPS);
            % bring in to right-canonical form
            testCase.testMPS = Can(testCase.testMPS, testCase.testMPSLength : -1 : 2, 'R');
            % build new state vector
            stateVec_R = Rebuild(testCase.testMPS);
            % assert they are the same up to some phase difference
            testCase.assertEqual(abs(stateVec_init), abs(stateVec_R), 'AbsTol', testCase.absTol);
        end
    end

end
