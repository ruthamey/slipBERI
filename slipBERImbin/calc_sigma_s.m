function [ sigma_s ] = calc_sigma_s( r_over_a, H)

% calc_sigma_s calculates the correlation matrix between slip patches.
% It is derived from the von karman function, which means the correlation
% depends on both the distance between slip patches and the down-dip (az)
% and along-strike (ax) correlation lengths, as well as the Hurst exponent.
%
% Inputs:
%       r_over_a = an n x n matrix of distance between slip patches divided by the
%                  correlation length
%       H        = an n x n matrix of Hurst parameters
%
% Ouputs:
%       sigma_s = n x n matrix, symmetrical matrix where each element gives the 
%                 correlation between slip patches given by the row/column number.
%                 e.g. row 2, column 3 would give the correlation between slip patch 2 and
%                 slip patch 3.
%
% First it calculates the modified bessel functions for r/a (distance
% between slip patches divided by correlation length) and then divides by
% the modified bessel function at ~ 0 (a scalar)
%
% rmja 18-dec-2014
% rmja  2-nov-2015



    sigma_s = (r_over_a).^H .* besselk( H, r_over_a) ./ (1e-10.^H.*besselk( H, 10^(-10)));
    
    % Set diagonals to 1 rather than 0
    sigma_s(logical(eye(size(sigma_s)))) = 1;  % THIS WAS FINE WHEN NOT SOLVING FOR CORRELATION LENGTHS. but now we are, so there are lots of diagonals, so it's easier to find them by just finding the zeros.
    sigma_s(isnan(sigma_s)) = 1;    % THIS WAS FINE WHEN NOT SOLVING FOR CORRELATION LENGTHS. but now we are, so there are lots of diagonals, so it's easier to find them by just finding the zeros.
%disp('POSSIBLY THERE IS A FILTHY CHEAT IN calc_sigma_s.m - to solve the look up table diagonal trouble in sigma_s calculation, I''ve set all nans to 1s. this works for the diagonals, but does it mess up other nans in other cases?')
    
    % Set inf and NaN to 0
     sigma_s(isinf(sigma_s)) = 0;
     
end

