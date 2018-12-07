function [] = calculate_moment(filename)
% calculate_moment returns moment calculates for a given saved Bayesian
% run.
% The code returns:
%       - Moment of the MAP solution
%       - Mean moment of all saved solutions
%       - 95% confidence intervals of the moments of all saved solutions
%
% R.M.J.Amey 2018

precision = 4;

load(filename, 'patch_MAP')
load(filename, 'elastic_params')
load(filename, 'spatial_model2column')
load(filename, 'spatial_model3column')

M0_MAP = sum((elastic_params.mu_okada) * spatial_model2column .* spatial_model3column .* (patch_MAP));
       
disp(['Moment of the MAP solution is ', num2str(M0_MAP, precision)])


% Calculate moment of ALL saved solutions, the mean of this, and 95% confidences
load(filename, 'slip_keep')

M0_all = sum((elastic_params.mu_okada) * spatial_model2column .* spatial_model3column .* (slip_keep));
mean_M0_all = mean(M0_all);
M0_conf_intervals = prctile(M0_all,[2.5 97.5]);

disp(['Mean M0 of all solutions is ', num2str(mean_M0_all, precision)])
disp(['95% confidence intervals are ', num2str(M0_conf_intervals(1), precision), ' and ', num2str(M0_conf_intervals(2), precision)])

figure;
hist(M0_all, 50);
xlabel('moment (Nm)')
ylabel('frequency')
hold on
plot([mean_M0_all mean_M0_all], [0 size(slip_keep,2)/50], 'r', 'LineWidth', 1)
plot([M0_MAP M0_MAP], [0 size(slip_keep,2)/50], 'g', 'LineWidth', 1)
plot([M0_conf_intervals(1) M0_conf_intervals(1)], [0 size(slip_keep,2)/50], 'y', 'LineWidth', 1)
plot([M0_conf_intervals(2) M0_conf_intervals(2)], [0 size(slip_keep,2)/50], 'y', 'LineWidth', 1)
legend('Hist', 'Mean', 'MAP', '95%', '95%')

end