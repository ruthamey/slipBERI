function [ H ] = calc_direction_dependent_H( H_dd, H_as, angle )

% calc_direction_dependent_H calculates the direction dependent Hurst
% parameters for each pair of n slip patches.
%
% Inputs:
%   H_dd = down-dip Hurst number (scalar)
%   H_as = along-strike Hurst number (scalar)
%   angle = n x n matrix of angles between slip patches
%
% Ouputs:
%   H = n x n matrix of direction dependent Hurst parameters
%
% For two patches along strike, the Hurst number will be H_as
% For two patches down dip, the Hurst number will be H_dd
% For patches which are some combination of distances along strike and
% down-dip of each other, this function gives a value between H_as and
% H_dd, as a function of the angle between them
%
% Ruth Amey 4-mar-2015

% Scale the angle, so it ranges from 0 to 1 ... recall that the atan function values from 0 to pi/2
scaled_angle = angle / (pi / 2);

% Using angle between, calculate direction dependent Hurst
H = H_dd *scaled_angle + H_as *(1-scaled_angle); % direction dependent Hurst

% Set diagonals to the average
H(logical(eye(size(H)))) = (H_dd + H_as) / 2;

end

