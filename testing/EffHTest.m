% EffHTest.m
% testing class for the effective Hamitlonian forming function EffH.m
% IMPORTANT: REQUIRES FILE CONTAINING THE MPO YOU WISH TO TEST AGAINST
% Oliver Thomson Brown
% 2015-08-03

classdef (SharedTestFixtures={matlab.unittest.fixtures.PathFixture('../dev/')}) EffHTest < matlab.unittest.TestCase

    properties
        absTol = 1E-14;
        testMPO = load('MPO_on-site+hop.mat', 'mpo');
        testHILBY = 2;
        testMPSLength = 9;
        testCOMPRESS = 0;
        testMPS;
        lBlock;
        rBlock;
    end

    methods (TestMethodSetup)
        function createMPS(testCase)
            testCase.testMPS = CompMPS(testCase.testHILBY, testCase.testMPSLength, testCase.testCOMPRESS);
        end
        
        function contractBlocks(testCase)
            testCase.lBlock = cell(testCase.testMPSLength, 1);
            testCase.lBlock{1} = ones(1,1,1);
            testCase.rBlock = cell(testCase.testMPSLength, 1);
            testCase.rBlock{testCase.testMPSLength} = ones(1,1,1);
            
            for site = 2 : 1 : testCase.testMPSLength
                testCase.lBlock{site} = GrowBlock(testCase.testMPS, testCase.testMPO.mpo, (site - 1), 'L', testCase.lBlock);
            end
            for site = (testCase.testMPSLength - 1) : -1 : 1
                testCase.rBlock{site} = GrowBlock(testCase.testMPS, testCase.testMPO.mpo, (site + 1), 'R', testCase.rBlock);
            end
        end
    end

    methods (Test)
        function hermiticityTest(testCase)
            % form effective Hamiltonian at each site and check it's Hermitian!
            for site = 1 : 1 : testCase.testMPSLength
                if site == 1
                    mpodex = 1;
                elseif site == testCase.testMPSLength
                    mpodex = 3;
                else
                    mpodex = 2;
                end  

                [rowSz, colSz, ~] = size(testCase.testMPS{site});
                effectiveHamiltonian = EffH(testCase.testHILBY, rowSz, colSz, testCase.lBlock{site}, testCase.testMPO.mpo{mpodex}, testCase.rBlock{site});

                % assert H - H' = 0
                effHSz = size(effectiveHamiltonian);
                hermTest = effectiveHamiltonian - ctranspose(effectiveHamiltonian);
                testCase.assertEqual(hermTest, zeros(effHSz), 'AbsTol', testCase.absTol);
            end
        end
    end
end
