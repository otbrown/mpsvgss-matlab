% ConvTest.m
% function which checks if the absolute value of the ground-state energy has been the same for the last $sampleSize site updates to within $threshold -- a simple convergence test
% Oliver Thomson Brown
% 2015-03-30
%
% [RETURN]
% convFlag	: test result, MATLAB doesn't use boolean values except in MuPad so convFlag == 1.0 if successful and convFlag == 0.0 if test is failed
%
% [INPUTS]
% energyTracker	: complex array, contains the energy of the system after each site update
% sampleSize	: double, the number of results to check for convergence
% threshold	: double, convergence threshold -- should be O(1E-14) or less

function [ convFlag ] = ConvTest(energyTracker, sampleSize, threshold)
	convFlag = 0;	

	if sampleSize >= length(energyTracker)
	else
		sample = energyTracker(end - sampleSize : end) - energyTracker(end);
		worstCase = max( abs(sample) );
		if worstCase < threshold
			convFlag = 1;
		end
	end			
end
