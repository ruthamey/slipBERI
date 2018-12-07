function [ prob ] = calc_logprior_laplacian( slip, L, alpha, n_fault_strands, first, last)
% 
% This code is really similar to calc_logprior_VK - see that for details
%
% Equation 3: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JB009017
% 
% Ruth Amey 5-jan-2016

prob = zeros(n_fault_strands,1);

    for i = 1 : n_fault_strands
        slip_temp = slip(first(i):last(i),1);
        L_temp = L(first(i):last(i),first(i):last(i));
        exponent = (-0.5 * (1/alpha(i)) * (L_temp*slip_temp)' * (L_temp* slip_temp));
        prob(i,1) = exponent;
    end

end