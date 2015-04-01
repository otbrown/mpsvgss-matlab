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
	convFlag = 0;
	updateCount = 0;
	direction = 'L';
	route = 1 : 1 : L;
	groundMPS = init_mps;
	groundMPS = Can(groundMPS, L : -1 : 2, 'R');

	% calculate initial state energy
	rightBlock = RBlock(groundMPS, mpo, 1);
	energyTracker = Expect(groundMPS, mpo, 1, rightBlock, 1);

	while updateCount < RUNMAX && ~convFlag
		for targetSite = route
			[rowMax, colMax, HILBY] = size(groundMPS{targetSite});

			if targetSite == 1
				mpodex = 1;
				groundMPS = Can(groundMPS, 2, 'R');
			elseif targetSite == L
				mpodex = 3;
				groundMPS = Can(groundMPS, L - 1, 'L');
			else
				mpodex = 2;
				groundMPS = Can(groundMPS, previousTarget, direction);
			end

			leftBlock = LBlock(groundMPS, mpo, targetSite);
			rightBlock = RBlock(groundMPS, mpo, targetSite);

			effectiveHamiltonian = EffH(HILBY, rowMax, colMax, leftBlock, mpo{mpodex}, rightBlock);

			[eigVec, energyTracker(end + 1)] = eigs(effectiveHamiltonian, 1, 'sr');

			groundMPS{targetSite} = SiteUpdate(eigVec, rowMax, colMax, HILBY);

			updateCount = updateCount + 1;		% updateCount++
			previousTarget = targetSite;

			fprintf('Update %u: E = %.5f\n', updateCount, real(energyTracker(end)));
			
			convFlag = ConvTest(energyTracker, L, THRESHOLD);

			% exit the loop if either the system is converged or the maximum number of updates has been reached
			if convFlag
				fprintf('Converged.\n');
				break;
			elseif updateCount >= RUNMAX
				fprintf('Failed to converge.\n');
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
