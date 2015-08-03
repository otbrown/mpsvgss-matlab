% MPSNormTest.m
% test for MPSNorm.m, a function which normalises a matrix product state
% Oliver Thomson Brown
% 2015-07-17

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev')})  MPSNormTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        testMPS;
    end

    properties (MethodSetupParameter)
        HILBY = {2, 2, 3, 4, 5, 2, 3, 4, 5};
        mpsLength = {5, 7, 6, 5, 5, 7, 6, 5, 5};
        COMPRESS = {0, 0, 0, 0, 0, 2, 6, 12, 15};
    end

    methods (TestMethodSetup, ParameterCombination='sequential')
        function MethodSetup(testCase, HILBY, mpsLength, COMPRESS)
            testCase.testMPS = CompMPS(HILBY, mpsLength, COMPRESS);
        end
    end

    methods (Test)
        function testNorm(testCase)
            psi = Rebuild(testCase.testMPS);
            norm = ctranspose(psi) * psi;
            testCase.assertEqual(norm, 1, 'AbsTol', testCase.absTol);
        end
    end
end         
