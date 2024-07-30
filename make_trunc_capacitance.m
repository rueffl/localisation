function [outputArg1,outputArg2] = make_trunc_capacitance(r,N,lij,alpha)
%MAKE_TRUNC_CAPACITANCE   Construct the truncated capacitance matrix
%   r:      truncation radius
%   N:      number of resonators in the unit cell
%   lij:    spacing between the resonators
%   alpha:  quasi periodicity

    Ir = 2*floor(r)+1; % number of lattice nodes with distance less than r to the origin

end