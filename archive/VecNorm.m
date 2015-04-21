% VecNorm.m
% script to normalise a vector
% ... I am lazy.
% Oliver Thomson Brown
% 20/10/2014

function [ vec ] = VecNorm(init)

N = sqrt( ctranspose(init) * init );
vec = init / N;

end
