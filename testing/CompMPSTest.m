% CompMPSTest.m
% tests for CompMPS.m function for generating arbitray matrix product states
% Oliver Thomson Brown
% 2015-07-17

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) CompMPSTest < matlab.unittest.TestCase

    properties
        testMPS;
        testMPSLength;
        testCOMPRESS;
        testHILBY;
    end

    properties (MethodSetupParameter)
        HILBY = {2, 2, 3, 4, 5, 2, 3, 4, 5};
        mpsLength = {5, 7, 6, 5, 5, 7, 6, 5, 5};
        COMPRESS = {0, 0, 0, 0, 0, 2, 6, 12, 15};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(testCase, HILBY, mpsLength, COMPRESS)
            testCase.testMPS = CompMPS(HILBY, mpsLength, COMPRESS);
            testCase.testMPSLength = mpsLength;
            testCase.testCOMPRESS = COMPRESS;
            testCase.testHILBY = HILBY;
        end
    end

    methods (Test)
        function testClass(testCase)
            testCase.assertClass(testCase.testMPS, 'cell');
        end

        function testCellShape(testCase)
            testCase.assertSize(testCase.testMPS, [testCase.testMPSLength, 1]);
        end

        function testComplex(testCase)
            for site = 1 : 1 : testCase.testMPSLength
                testCase.assertFalse(isreal(testCase.testMPS{site}));
            end
        end

        function testMatShape(testCase)
            if testCase.testCOMPRESS == 0
                testCase.testCOMPRESS = Inf;
            end
            % initialise rowSz and colSz
            rowSz = 1;
            colSz = min(testCase.testHILBY, testCase.testCOMPRESS);
            % sweep through each site in the matrix product state
            for site = 1 : 1 : testCase.testMPSLength;
                testCase.assertSize(testCase.testMPS{site}, [rowSz, colSz, testCase.testHILBY]);
                % gauge whether there are more sites to the left, or to the right of the current one
                lLen = site + 1;
                rLen = testCase.testMPSLength - site - 1;
                len = min(lLen, rLen);
                % set next site rowSz and colSz
                rowSz = colSz;
                colSz = min(testCase.testCOMPRESS, testCase.testHILBY^len);
            end
        end
    end
end
