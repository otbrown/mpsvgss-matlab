% CanTest.m
% a script to test the function Can.m or more specifically,
% its sub-functions LCan.m and RCan.m
% Oliver Thomson Brown
% 15-07-20

% invoke empty interactive test-case object, to allow assertEqual 
tc = matlab.unittest.TestCase.forInteractiveUse;

% set absolute tolerance
absTol = 1E-14;

% test mps
mps250 = CompMPS(2, 5, 0);
mps270 = CompMPS(2, 7, 0);
mps360 = CompMPS(3, 6, 0);
mps450 = CompMPS(4, 5, 0);
mps550 = CompMPS(5, 5, 0);
mps272 = CompMPS(2, 7, 2);
mps366 = CompMPS(3, 6, 6);
mps4512 = CompMPS(4, 5, 12);
mps5515 = CompMPS(5, 5, 15);

% left-normalise mps
mps250 = Can(mps250, 1 : 4, 'L');
mps270 = Can(mps270, 1 : 6, 'L');
mps360 = Can(mps360, 1 : 5, 'L');
mps450 = Can(mps450, 1 : 4, 'L');
mps550 = Can(mps550, 1 : 4, 'L');
mps272 = Can(mps272, 1 : 6, 'L');
mps366 = Can(mps366, 1 : 5, 'L');
mps4512 = Can(mps4512, 1 : 4, 'L');
mps5515 = Can(mps5515, 1 : 4, 'L');

classdef CanTest < matlab.unittest.TestCase
