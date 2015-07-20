% ConvTestTest.m
% somewhat awkwardly named test for ConvTest.m which determines whether or not the ground-state search has converged
% Oliver Thomson Brown
% 2015-07-17

% test energyTracker arrays
eT20 = ones(20,1);
eT20m12 = ones(20,1) + 1E-12 * rand(20,1);

%% assert-true
assert(ConvTest(eT20,5,0.1) == 1, 'Error: ConvTest reports a false negative. Code: eT205.1');
assert(ConvTest(eT20,10,0.1) == 1, 'Error: ConvTest reports a false negative. Code: eT2010.1');
assert(ConvTest(eT20,19,0.1) == 1, 'Error: ConvTest reports a false negative. Code: eT2019.1');
assert(ConvTest(eT20m12,5,1E-12) == 1, 'Error: ConvTest reports a false negative. Code: eT20m125m12');
assert(ConvTest(eT20m12,10,1E-12) == 1, 'Error: ConvTest reports a false negative. Code: eT20m1210m12');
assert(ConvTest(eT20m12,19,1E-12) == 1, 'Error: ConvTest reports a false negative. Code: eT20m1219m12');

%% assert-false
assert(ConvTest(eT20,25,0.1) == 0, 'Error: ConvTest reports a false positive on sample size overflow. Code: eT2025.1');
assert(ConvTest(eT20,30,0.1) == 0, 'Error: ConvTest reports a false positive on sample size overflow. Code: eT2030.1');
assert(ConvTest(eT20,100,0.1) == 0, 'Error: ConvTest reports a false positive on sample size overflow. Code: eT20100.1');
assert(ConvTest(eT20m12,5,1E-13) == 0, 'Error: ConvTest reports a false positive above threshold. Code: eT20m125m13');
assert(ConvTest(eT20m12,10,1E-13) == 0, 'Error: ConvTest reports a false positive above threshold. Code: eT20m1210m13');
assert(ConvTest(eT20m12,19,1E-13) == 0, 'Error: ConvTest reports a false positive above threshold. Code: eT20m1219m13');
assert(ConvTest(eT20m12,5,1E-16) == 0, 'Error: ConvTest reports a false positive above threshold. Code: eT20m125m16');
assert(ConvTest(eT20m12,10,1E-16) == 0, 'Error: ConvTest reports a false positive above threshold. Code: eT20m1210m16');
assert(ConvTest(eT20m12,19,1E-16) == 0, 'Error: ConvTest reports a false positive above threshold. Code: eT20m1219m16');
