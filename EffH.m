% EffH.m
% Oliver Thomson Brown
% 2015-03-02
% DOCSTRING!

%{
close all
clear all
clc

HILBY = 2;
L = 5;
COMPRESS = 0;
TARGET = 1;

mps = CompMPS(HILBY, L, COMPRESS);

% GENERATE AN MPO

om = 20;
J = 5;

I2 = eye(2);                    % basis blocks
Z2 = zeros(2);
spinUp = [0, 0; 1, 0];
spinDown = [0, 1; 0, 0];
spinCount = spinUp * spinDown;
spinZi = spinUp * spinDown - spinDown * spinUp;

W1 = [ om * spinCount, -J * spinDown, -J * spinUp, I2 ];        % mpo matrices, constructed from block matrices
Wi = [ I2, Z2, Z2, Z2; spinUp, Z2, Z2, Z2; spinDown, Z2, Z2, Z2; om * spinCount, -J * spinDown, -J * spinUp, I2 ];
WL = [ I2; spinUp; spinDown; om * spinCount ];

%W1 = [ om * spinCount, J * spinZi, I2 ];
%Wi = [ I2, Z2, Z2; spinZi, Z2, Z2; om * spinCount, J * spinZi, I2];
%WL = [ I2; spinZi; om * spinCount];

mpo = cell(3,1);

mpo{1} = W1;
mpo{2} = Wi;
mpo{3} = WL;

impo = cell(3,1);               % Identity MPO
impo{1} = I2;
impo{2} = I2;
impo{3} = I2;

% set up for effh

leftBlock = LBlock(mps, mpo, TARGET);
rightBlock = RBlock(mps, mpo, TARGET);

[rowMax, colMax, HILBY] = size( mps{TARGET} );

if TARGET == L
	mpodex = 3;
elseif TARGET == 1
	mpodex = 1;
else
	mpodex = 2;
end

mpo = mpo{mpodex};
%}

function [ effectiveHamiltonian ] = EffH(HILBY, rowMax, colMax, leftBlock, mpo, rightBlock);
	% gather data	
	opRowMax = ( size(mpo, 1) / HILBY ) - 1;
	opColMax = ( size(mpo, 2) / HILBY ) - 1;

	% pre-allocate
	dimension = HILBY * rowMax * colMax;
	effectiveHamiltonian = sparse(dimension, dimension);

	% LOOP THE LOOP
	for braState = 0 : 1 : HILBY - 1
		for ketState = 0 : 1 : HILBY - 1
			for row = 0 : 1 : rowMax - 1
				for col = 1 : 1 : colMax
					for conjRow = 1 : 1 : colMax 
						for conjCol = 0 : 1 : rowMax - 1
							% joint indexing
							jRow = braState * colMax * rowMax + conjCol * colMax + conjRow;
							jCol = ketState * rowMax * colMax + row * colMax + col;
							% introduce block indices and contract LWR
							for opRow = 0 : 1 : opRowMax
								for opCol = 0 : 1 : opColMax
									effectiveHamiltonian(jRow, jCol) = effectiveHamiltonian(jRow, jCol) + ...
														leftBlock(conjCol + 1, row + 1, opRow + 1) ...
														* mpo(opRow * HILBY + braState  + 1, opCol * HILBY + ketState + 1) ...
														* rightBlock(conjRow, col, opCol + 1);
									% END TIMES
								end
							end
						end
					end
				end
			end
		end
	end
	 
end
