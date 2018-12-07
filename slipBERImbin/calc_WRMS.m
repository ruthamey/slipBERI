function [ WRMS ] = calc_WRMS( slip, G, d, inv_sigma_d, offset )

% calc_rms calculates the residuals between the observed data and predicted
% data for the given slip values.
%
% Inputs: (slip, G, d, inv_sigma_d, offset)
%
% This function calculates the forward model, by multiplying the given slip
% vector by the appropriate matrix of Green's functions. This gives a
% forward model, or predicted set of observations. The difference between
% this forward model and the observed observations is calculated.
%
% Ruth Amey 24-feb-2015
% Ruth Amey 22-oct-2015  Updated to WRMS not RMS

WRMS = ((d - G*slip)+offset)' * inv_sigma_d * ((d - G*slip)+offset);

end

