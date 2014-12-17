% Can.m
% accessor function for the left- and right-canonisation functions,
% LCan and RCan; in particular Can ensures the supplied route is valid
% Oliver Thomson Brown
% 13/11/2014
%
% [RETURN]
% cmps		: cell array containing canonised matrix product state, A_state_site = cmps{site}(:,:,state)
%
% [INPUT]
% mps		: cell array containing a matrix product state, A_state_site = cmps{site}(:,:,state)
% route		: vector which contains sites which are to be canonised
% DIRECTION	: the direction in which the state should be canonised; 'L' or 'R'

function [cmps] = Can(mps, route, DIRECTION)
	
	if DIRECTION == 'L'
		L = size(mps,1);
		if route(end) == L
			route = route(1 : end - 1);
		end
		cmps = LCan(mps, route);
	elseif DIRECTION == 'R'
		if route(end) == 1
			route = route(1 : end - 1);
		end
		cmps = RCan(mps, route);
	else
		fprintf('Please enter a valid canonisation direction, ''L'' or ''R''.\n');
	end
end
