function [ exponent ] = calc_loglikely( slip, d, G, inv_sigma_d, offset, beta )
% For a given slip distribution, this function looks at the residuals
% between observed surface displacement and those predicted for the
% given slip distribution (calculated by G*s), where G is the kernel of
% surface displacements calculated from Okada 1985
%
% Inputs:
%     slip    = slip on each fault patch
%      d      = data points (InSAR, GPS, etc)
%      G      = correct kernel, that relates slip on each patch to each datapoint
% inv_sigma_d = inverse of variance-covariance matrix of data and model errors
%     offset  = InSAR offset (values for each InSAR datapoint, 0 for any GPS datapoints)
%     beta    = data hyperparameter
%
% Output:
%   exponent = exponent of the likelihood function
%
% NOTE: this calculates the log of the likelihood function, since the
% likelihood function is otherwise too small.
% The full likelihood equation is:
% likelihood = (2*pi)^(N/2) * det_sigma_d^(1/2) * exp ( (-0.5 * (1/beta^2) * ((d-G*slip)+offset)' * inv_sigma_d * ((d-G*slip)+offset) ))
% But here we only calculate the exponent. Thus when considering ratios in
% the main script, instead of calculating the ratio as:
%           likelihood ratio = prior_trial / prior_curr;
% We calculate it as:
%       exponent_ratio = exp( logprior_trial - logprior_curr);
%
% Ruth Amey 11-dec-2014   Written
% Ruth Amey 26-aug-2015   Better explained
% Ruth Amey  7-oct-2015   Add nan check

if length(beta) == 0
   beta = 1; 
end

if length(offset) == 0
   offset = zeros( length(d),1); 
end

exponent = (-0.5 * (1/beta^2) * (d-(G*slip+offset))' * inv_sigma_d * (d-(G*slip+offset)));

if isnan(exponent)
    disp('DANGER: A rogue NaN has been located in calc_loglikely')
    keyboard
end

end

