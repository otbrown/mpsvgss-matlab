% EffH.m
% Oliver Thomson Brown
% 2015-03-02
% DOCSTRING!

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
					for conjRow = 0 : 1 : colMax - 1
						for conjCol = 1 : 1 : rowMax
							% joint indexing
							jRow = braState * colMax * rowMax + conjRow * rowMax + conjCol;
							jCol = ketState * rowMax * colMax + row * colMax + col;
							% introduce block indices and contract LWR
							for opRow = 0 : 1 : opRowMax
								for opCol = 0 : 1 : opColMax
									effectiveHamiltonian(jRow, jCol) = effectiveHamiltonian(jRow, jCol) + ...
														leftBlock(conjRow + 1, row + 1, opRow + 1) ...
														* mpo(opRow * HILBY + braState  + 1, opCol * HILBY + ketState + 1) ...
														* rightBlock(conjCol, col, opCol + 1);
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
