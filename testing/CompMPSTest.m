% CompMPSTest.m
% tests for CompMPS.m function for generating arbitray matrix product states
% Oliver Thomson Brown
% 2015-07-17

% test matrix product states
mps270 = CompMPS(2, 7, 0);      % uncompressed
mps360 = CompMPS(3, 6, 0);
mps450 = CompMPS(4, 5, 0);
mps550 = CompMPS(5, 5, 0);
mps272 = CompMPS(2, 7, 2);      % compressed
mps366 = CompMPS(3, 6, 6);
mps4512 = CompMPS(4, 5, 12);
mps5515 = CompMPS(5, 5, 15);

% absolute tolerance 
% no actual floating point testing here because we're just checking size, shape, and types, so tol = 0
tol = 0;

% test outputs are correct type
assert(isa(mps270,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');
assert(isa(mps360,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');
assert(isa(mps450,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');
assert(isa(mps550,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');
assert(isa(mps272,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');
assert(isa(mps366,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');
assert(isa(mps4512,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');
assert(isa(mps5515,'cell'), 'Fundamental problem: CompMPS is not returning cell arrays.');

% test outputs are correct shape
sz270 = size(mps270);
assert(sz270(1) == 7, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz270(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
sz360 = size(mps360);
assert(sz360(1) == 6, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz360(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
sz450 = size(mps450);
assert(sz450(1) == 5, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz450(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
sz550 = size(mps550);
assert(sz550(1) == 5, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz550(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
sz272 = size(mps272);
assert(sz272(1) == 7, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz272(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
sz366 = size(mps366);
assert(sz366(1) == 6, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz366(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
sz4512 = size(mps4512);
assert(sz4512(1) == 5, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz4512(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
sz5515 = size(mps5515);
assert(sz5515(1) == 5, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');
assert(sz5515(2) == 1, 'Fundamental problem: CompMPS is returning cell arrays of the wrong shape.');

% test CompMPS is providing a complex-valued MPS
for site = 1 : 1 : length(mps270)
    A = mps270{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end
for site = 1 : 1 : length(mps360)
    A = mps360{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end
for site = 1 : 1 : length(mps450)
    A = mps450{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end
for site = 1 : 1 : length(mps550)
    A = mps550{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end
for site = 1 : 1 : length(mps272)
    A = mps272{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end
for site = 1 : 1 : length(mps366)
    A = mps366{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end
for site = 1 : 1 : length(mps4512)
    A = mps4512{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end
for site = 1 : 1 : length(mps5515)
    A = mps5515{site};
    assert(~isreal(A), 'Error: CompMPS is returning real-valued matrices');
end

% test the matrices in the mps are consistently sized -- the method used here is complex, but definitely correct
testSizes(mps270, 2, 0);
testSizes(mps360, 3, 0);
testSizes(mps450, 4, 0);
testSizes(mps550, 5, 0);
testSizes(mps272, 2, 2);
testSizes(mps366, 3, 6);
testSizes(mps4512, 4, 12);
testSizes(mps5515, 5, 15);
