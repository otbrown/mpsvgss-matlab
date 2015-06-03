% Ground.m
% function which performs a ground-state search given an initial Matrix Product State and a Hamiltonian in mpo form
% Oliver Thomson Brown
% 2015-03-31
%
% [RETURN]
% groundMPS	: cell array, the ground state of the system in MPS form
% energyTracker : contains the eigenvalue from each site update, ground state energy is energyTracker(end)
%
% [INPUTS]
% init_mps	: cell array, initial matrix product state -- must be the right size and shape, but otherwise arbitrary **ASSUMED TO BE NORMALISED!**
% mpo		: 3 x 1 cell array, the Hamiltonian of the system in matrix product operator form
% THRESHOLD	: double, convergence threshold
% RUNMAX	: double (integer), the maximum number of updates the code should perform before giving up

function [ groundMPS, energyTracker ] = Ground(init_mps, mpo, THRESHOLD, RUNMAX)
	% gather data from input args
	L = length(init_mps);

	% initialise variables
	firstConvFlag= 0;
	fullConvFlag = 0;
	updateCount = 0;
	direction = 'L';
	route = 1 : 1 : L;
	groundMPS = init_mps;
	groundMPS = Can(groundMPS, L : -1 : 2, 'R');

	% build left and right contractions
	left = cell(L, 1);
	right = cell(L, 1);
	left{1} = ones(1,1,1);
	right{L} = ones(1,1,1);
	
	fprintf('Building contraction from right-hand side.\n');
	for targetSite = L - 1 : -1 : 1
		right{targetSite} = GrowBlock(groundMPS, mpo, targetSite + 1, 'R', right);
	end
	fprintf('Built.\n');

	% calculate initial state energy
	energyTracker = Expect(groundMPS, mpo, 1, right{1}, 1);

	while updateCount < RUNMAX && ~fullConvFlag
		for targetSite = route
			[rowMax, colMax, HILBY] = size(groundMPS{targetSite});

			if targetSite == 1
				mpodex = 1;
				groundMPS = Can(groundMPS, 2, 'R');
				right{1} = GrowBlock(groundMPS, mpo, 2, 'R', right);
			elseif targetSite == L
				mpodex = 3;
				groundMPS = Can(groundMPS, L - 1, 'L');
				left{L} = GrowBlock(groundMPS, mpo, L - 1, 'L', left);
			else
				mpodex = 2;
				groundMPS = Can(groundMPS, previousTarget, direction);
				if direction == 'L'
					left{targetSite} = GrowBlock(groundMPS, mpo, previousTarget, direction, left);
				elseif direction == 'R'
					right{targetSite} = GrowBlock(groundMPS, mpo, previousTarget, direction, right);
				end
			end

			leftBlock = left{targetSite};
			rightBlock = right{targetSite};

			effectiveHamiltonian = EffH(HILBY, rowMax, colMax, leftBlock, mpo{mpodex}, rightBlock);
			
			[eigVec, energyTracker(end + 1)] = eigs(effectiveHamiltonian, 1, 'sr');

			groundMPS{targetSite} = SiteUpdate(eigVec, rowMax, colMax, HILBY);

			updateCount = updateCount + 1;		% updateCount++
			previousTarget = targetSite;

			fprintf('Target Site: %u\n', targetSite);
			fprintf('Update %u: E = %.5f\n', updateCount, real(energyTracker(end)));
			if ~firstConvFlag
				firstConvFlag = ConvTest(energyTracker, 5, 1E-6);
				if firstConvFlag
					filename = sprintf('%dx%dGSSchkpnt', L, HILBY);
					save(filename, 'init_mps', 'groundMPS', 'energyTracker', '-v7.3');
					fprintf('First convergence threshold reached, checkpoint created.\n');
				end
			else
				fullConvFlag = ConvTest(energyTracker, L, THRESHOLD);
			end

			% exit the loop if either the system is converged or the maximum number of updates has been reached
			if fullConvFlag
				fprintf('Converged.\n');
				break;
			elseif updateCount >= RUNMAX
				fprintf('Failed to converge down to %d.\n', THRESHOLD);
				break;
			end
		end

		% flip it and reverse it
		if direction == 'L'
			direction = 'R';
		elseif direction == 'R'
			direction = 'L';
		end

		route = flip(route);
	end
end
