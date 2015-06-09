% ParGround.m
% function which performs a ground-state search given an initial Matrix Product State and a Hamiltonian in mpo form
% Oliver Thomson Brown
% 2015-06-24
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
% NUM_WORKERS : double (integer), the number of workers to assign

function [ groundMPS, energyTracker ] = ParGround(init_mps, mpo, THRESHOLD, RUNMAX, NUM_WORKERS)
tStart = tic;

% gather data from input args
L = length(init_mps);

% initialise variables
firstConvFlag = 0;
secondConvFlag = 0;
fullConvFlag = 0;
updateCount = 0;
groundMPS = init_mps;

% build left and right contractions
left = cell(L, 1);
right = cell(L, 1);
left{1} = ones(1,1,1);
right{L} = ones(1,1,1);

groundMPS = Can(groundMPS, L, 'R');

fprintf('Building right and left contraction blocks.\n');
% taking advantage of fact that MPSNorm left-normalises
for targetSite = 2 : 1 : L
    left{targetSite} = GrowBlock(groundMPS, mpo, targetSite - 1, 'L', left);
end

groundMPS = Can(groundMPS, L:-1:2, 'R');
for targetSite = L - 1 : -1 : 1;
    right{targetSite} = GrowBlock(groundMPS, mpo, targetSite + 1, 'R', right);
end
fprintf('Built.\n');

% calculate initial state energy
energyTracker = Expect(groundMPS, mpo, 1, right{1}, 1);

parpool(NUM_WORKERS);

while updateCount < RUNMAX && ~fullConvFlag
    parfor targetSite = 1 : 1 : L
        [rowMax, colMax, HILBY] = size(groundMPS{targetSite});
        
        if targetSite == 1
            mpodex = 1;
        elseif targetSite == L
            mpodex = 3;
        else
            mpodex = 2;
        end
        
        leftBlock = left{targetSite};
        rightBlock = right{targetSite};
        
        effectiveHamiltonian = EffH(HILBY, rowMax, colMax, leftBlock, mpo{mpodex}, rightBlock);
                       
        [eigVec, ~] = eigs(effectiveHamiltonian, 1, 'sr');
        
        groundMPS{targetSite} = SiteUpdate(eigVec, rowMax, colMax, HILBY);
    end
    
    fprintf('Chain updated.\nRenormalising.\n')
    updateCount = updateCount + 1;
    
    % renormalisation, and reconstruction of left and right blocks
    groundMPS = MPSNorm(groundMPS);
    fprintf('Building right and left contraction blocks.\n');
    % taking advantage of fact that MPSNorm left-normalises
    for targetSite = 2 : 1 : L
        left{targetSite} = GrowBlock(groundMPS, mpo, targetSite - 1, 'L', left);
    end
    groundMPS = Can(groundMPS, L:-1:2, 'R');
    for targetSite = L - 1 : -1 : 1
        right{targetSite} = GrowBlock(groundMPS, mpo, targetSite + 1, 'R', right);
    end
    fprintf('Built.\n');

    energyTracker(end+1) = Expect(groundMPS, mpo, 1, right{1}, 1);
    
    fprintf('Update %u: E = %d\n', updateCount, energyTracker(end));
    
    HILBY = size(groundMPS, 3);

    if ~firstConvFlag
        firstConvFlag = ConvTest(energyTracker, 3, 1E-5);
        if firstConvFlag
            filename = sprintf('%dx%dGSSchkpnt', 3, HILBY);
            save(filename, 'init_mps', 'groundMPS', 'energyTracker', '-v7.3');
            fprintf('Convergence to 1E-5 reached, checkpoint created.\n');
            toc(tStart);
        end
    elseif ~secondConvFlag
        secondConvFlag = ConvTest(energyTracker, 3, 1E-7);
        if secondConvFlag
            fprintf('Convergence to 1E-7 reached.\n');
            toc(tStart);
        end
    else
        fullConvFlag = ConvTest(energyTracker, 3, THRESHOLD);
    end
    
    % exit the loop if either the system is converged or the maximum number of updates has been reached
    if fullConvFlag
        fprintf('Converged.\n');
        toc(tStart);
        break;
    elseif updateCount >= RUNMAX
        fprintf('Failed to converge down to %d.\n', THRESHOLD);
        toc(tStart);
        break;
    end
    
end

delete(gcp('nocreate'));

end


