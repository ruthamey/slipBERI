function [] = judge_convergence(filename)
%% JUDGEMENT DAY IS UPON US
%
%


fontsize_plot = 32;
set(0,'DefaultAxesFontSize',fontsize_plot)
set(0,'defaultAxesFontName', 'Times New Roman')

%% Histograms

%% 2D histograms

%% posterior

 load(filename, 'logposterior_keep')
 figure
 plot( 1:length(logposterior_keep), logposterior_keep)
 xlabel('Iterations')
 ylabel('Unscaled Probability')
 
 %% moment
 %can you not plot evolution of 95% confidence for some relevant parameter, like moment?
 
load(filename, 'patch_MAP')
load(filename, 'elastic_params')
load(filename, 'spatial_model2column')
load(filename, 'spatial_model3column')
load(filename, 'slip_keep')

M0_all = sum((elastic_params.mu_okada) * spatial_model2column .* spatial_model3column .* (slip_keep));

intervals_to_calc_95_conf_at = 1000:1000:length(slip_keep);
pcent95 = zeros(length(intervals_to_calc_95_conf_at),2);

%CLUNKY
for i = 1: length(intervals_to_calc_95_conf_at)
   pcent95(i,:) = prctile(M0_all( 1:intervals_to_calc_95_conf_at(i))',[2.5 97.5]);
end

figure;
plot(1:size(pcent95,1), pcent95(:,1), 'g', 'linewidth', 5)
hold on
plot(1:size(pcent95,1), pcent95(:,2), 'r', 'linewidth', 5)
xlabel('Iterations')
ylabel('Moment (Nm)')
 
 %% results of first third vs. results of last third

end




