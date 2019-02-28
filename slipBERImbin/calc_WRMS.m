function [ WRMS ] = calc_WRMS( slip, G, d, inv_sigma_d, offset )

% calc_WRMS calculates the residuals between the observed data and predicted
% data for the given slip values.
%
% Inputs: (slip, G, d, inv_sigma_d, offset)
%
% This function calculates the forward model, by multiplying the given slip
% vector by the appropriate matrix of Green's functions. This gives a
% forward model, or predicted set of observations. The difference between
% this forward model and the observed observations is calculated, weighted
% by the variance-covariance, normalised by the length and then squared.
% The square root is then taken to give thinverte WRMS.
%
% Ruth Amey 24-feb-2015
% Ruth Amey 22-oct-2015  Updated to WRMS not RMS

  residuals = (d - (G*slip+offset));
%  weighted_residuals_sq = residuals'* inv_sigma_d * residuals;
%  sum_weighted_residuals_sq =  sum(weighted_residuals_sq);
%  normalised_sum_weighted_residuals_sq = sum_weighted_residuals_sq / length(d);
%  sqrt_sum_weighted_residuals_sq = sqrt(normalised_sum_weighted_residuals_sq);

WRMS = sqrt(         (  sum(  residuals' * inv_sigma_d * residuals   )   /   length(d)  )         );

end

