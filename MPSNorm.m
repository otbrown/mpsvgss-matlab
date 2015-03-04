% MPSNorm.m
% Function which returns the (left-)normalised matrix product state.
% Oliver Thomson Brown
% 2015-03-04
% DOCSTRING!

function [ normalMPS ] = MPSNorm( mps )
	threshold = 1E-14;

	L = size( mps, 1 );			% TASHA YAR
	normalMPS = mps;

	for site = 1 : 1 : L - 1
		M = cat(1, normalMPS{site}(:, :, 1), normalMPS{site}(:, :, 2));
		[U, S, V] = svd(M);

		U(abs(U) < threshold) = 0;	% on the regz
		S(abs(S) < threshold) = 0;
		V(abs(V) < threshold) = 0;

		Urow = size(U, 1);
		rowLim = size(M, 1) / 2;
		colLim = size(M, 2);

		normalMPS{site}(:, :, :) = cat(3, U(1 : rowLim, 1 : colLim), U(Urow/2 + 1 : Urow/2 + rowLim, 1 : colLim));

		chain = S * ctranspose(V);	% listen to the wind blow

		rowLim = colLim;
		colLim = size(normalMPS{site + 1}(:, :, 1), 2);

		N = cat(3, chain * normalMPS{site + 1}(:, :, 1), chain * normalMPS{site + 1}(:, :, 2));

		normalMPS{site + 1} = cat(3, N(1 : rowLim, 1 : colLim, 1), N(1 : rowLim, 1 : colLim, 2));
	end

	M = cat(1, normalMPS{L}(:, :, 1), normalMPS{L}(:, :, 2));
	[U, S, V] = svd(M);

	U(abs(U) < threshold) = 0;

	Urow = size(U, 1);
	rowLim = size(M, 1) / 2;
	colLim = size(M, 2);

	normalMPS{L}(:, :, :) = cat(3, U(1 : rowLim, 1 : colLim), U(Urow/2 + 1 : Urow/2 + rowLim, 1 : colLim));		
end 
