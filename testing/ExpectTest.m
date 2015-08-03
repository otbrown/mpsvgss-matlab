% ExpectTest.m
% class testing for Expect.m, checks that the expectation of various MPS with themselves is 1
% for that reason it somewhat relies on the normalisation being correct
% Oliver Thomson Brown
% 2015-08-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')}) ExpectTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14
        iMPO;
        testMPS;
        testMPSLength;
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
            testCase.iMPO = [{eye(HILBY)}; {eye(HILBY)}; {eye(HILBY)}];
        end
    end

    methods (Test)
        function normExpectTest(testCase)
            % set target site about which expectation will be taken
            TARGET = ceil(testCase.testMPSLength / 2);
            % allocate contractions up to the target site
            lBlock = cell(testCase.testMPSLength);
            rBlock = cell(testCase.testMPSLength);
            lBlock{1} = ones(1,1,1);
            rBlock{testCase.testMPSLength} = ones(1,1,1);
            % compute contractions up to the target site
            for site = 2 : 1 : TARGET
                lBlock{site} = GrowBlock(testCase.testMPS, testCase.iMPO, (site - 1), 'L', lBlock);
            end
            for site = (testCase.testMPSLength - 1) : -1 : TARGET
                rBlock{site} = GrowBlock(testCase.testMPS, testCase.iMPO, (site + 1), 'R', rBlock);
            end
            % compute expectation value -- should be the norm of the state (which should be one!)
            expectationValue = Expect(testCase.testMPS, testCase.iMPO, lBlock{TARGET}, rBlock{TARGET}, TARGET);
            % assert
            testCase.assertEqual(expectationValue, 1, 'AbsTol', testCase.absTol);
        end
    end
end
