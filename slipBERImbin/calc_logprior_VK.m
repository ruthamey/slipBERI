function [ exponent, singularflag ] = calc_logprior_VK( slip, inv_sigma_s, alpha, n_fault_strands_for_smoothin, first, last)
%
% calc_logprior_VK calculates the VK prior probability of a slip
% distribution
%
% The calculation is weighted by the inverse of the correlation matrix 
% (inv_sigma_s), which is derived from the von karman function. It is also 
% weighted by the prior hyperparameter, alpha.
%
% We calculate the prior probability separately for each separately smoothed
% fault strands, so for n fault strands, the output have n rows
%
% Input:
%                   slip               = slip in each patch
%                   inv_sigma_s        = inverse of VK correlation matrix
%                   alpha              = prior hyperparameter
%        n_fault_strands_for_smoothing = the number of separately smoothed fault strands
%                   first              = first patch on each fault strand
%                   last               = last patch on each fault strand
%
% Output:
%               exponent term          =   (-0.5 * (1/alpha(i)) * s.' * inv_sigma_s_temp * s)   - one value per fault strand.
%
%
% NOTE: This calculates the logs of the prior likelihoods.
%
% Ruth Amey 27-feb-2014

exponent = zeros(n_fault_strands_for_smoothin,1);

for i = 1 : n_fault_strands_for_smoothin                                  % Calculate one term per smoothed fault strand
       
    slip_temp = slip(first(i):last(i),1);                                 % Select the right slip patches for that fault strand.
    inv_sigma_s_temp = inv_sigma_s(first(i):last(i),first(i):last(i));    % Select the right inv_sigma_s matrix for that fault strand.
    mean_slip = mean(slip_temp,1);
    s = slip_temp - mean_slip;
    exponent(i,1) = (-0.5 * (1/alpha(i)) * s.' * inv_sigma_s_temp * s);
    inv_sigma_s_temp = [];

end

if isnan(exponent)
    disp('DANGER: nan in calc_logprior_VK');
    keyboard
end

if exponent>0
    disp('DANGER: exponent in calc_logprior_VK is positive');
    %keyboard
    %disp('setting VK prior probability to 0');
    exponent = -1000000;
    singularflag = 1;
else
    singularflag = 0;
end

end