%  Script to plot the slip result from slipBERI.
%
% INPUTS, 'display' is a structure detailing what you'd like to display:
%         plotmean = String, 'yes' or 'no'. whether you'd like to plot the mean and standard deviation or not.
%         plotmode = String, 'yes' or 'no'. whether you'd like to plot the mode and standard deviation or not.
%         plotmedian = String, 'yes' or 'no'. whether you'd like to plot the median and standard deviation or not.
%         plotallsips = String, 'yes' or 'no'. whether you'd like one plot with the mean, mode, median, max likelihood (and true, if if in testing mode).
%         plotprob = String, 'yes' or 'no'. whether you'd like to plot how the probability changes with the number of iterations.
%         plothists = String. how many histograms you'd like to display. 'plothistall' means plot histograms for all the slip patches, 'plothistsample' means select some randomly to display
%         plotsurfacedisp = String, 'yes' or 'no'. Display surface displacement - arrows for GPS, colour for InSAR LOS.
%         plotmarginalPDFs = String, 'yes' or 'no'. this chooses the six patches with the highest slip and plots the marginal PDFs for this.
%         plot_resolution_matrix = String. 'yes' or 'no'.
%         plot3d = String. 'yes' or 'no' for whether to plot the fault geometry with slip mode.
%         plotMAP = String. 'yes' or 'no' for whether to plot the MAP (max a posteriori) and standard deviation or not
%         calc_confidence = String. 'yes' or 'no' for whether you want to plot the 95% confidence intervals for each patch.
%
% R.M.J.Amey 2018



%% Load data, if you need to

fontsize_plot = 32;
set(0,'DefaultAxesFontSize',fontsize_plot)
set(0,'defaultAxesFontName', 'Times New Roman')

oldversion = exist('total_n_slip_patches_true');
        if oldversion  == 0
            total_n_slip_patches = total_n_slip_patches;
        end
        
%         oldversion = exist('n_down_dip_patches_true');
%         if oldversion  == 0
%             total_n_down_dip_patches = n_down_dip_patches_for_smoothing;
%         end
        
        oldversion = exist('n_along_strike_patches_true');
        if oldversion  == 0
            n_along_strike_patches = n_along_strike_patches_for_smoothing;
        end

%% Find mean, mode, everything

if strcmp(invert.inversion_type, 'bayesian') == 1
    
    onlyonkeptpatches = slip_keep~=0;
    
        % Find mean slip on each slip patch
%         if strcmp(invert.solve_for_fault_size, 'yes') == 1
%             patch_mean = zeros(n_slip_patches_on_each_fault_strand,1);
%             rake_mean = zeros(n_slip_patches_on_each_fault_strand,1);
%             for i = 1:total_n_slip_patches
%                 patch_mean(i,1) = mean(slip_keep(i,onlyonkeptpatches(i,:)));
%                 rake_mean(i,1) = mean(rake_keep(i,onlyonkeptpatches(i,:)));
%             end
%             patch_mean(isnan(patch_mean)==1) = 0;
%             rake_mean(isnan(rake_mean)==1) = 0;
%         else
            patch_mean = mean(slip_keep, 2);
            rake_mean = mean(rake_keep, 2);
        %end
        %total_n_down_dip_patches = n_down_dip_patches(1);
        %total_n_along_strike_patches = sum(n_along_strike_patches);
        oldversion = exist('variable_patches_along_strike_with_depth');
        if oldversion  == 0
            variable_patches_along_strike_with_depth = 'no';
        end


        oldversions = exist('onoffidentifyerstoretotal');
        if oldversion  == 0
            onoffidentifyerstoretotal = onoffidentifyerkeep;
        end

        oldversion = exist('data.atolls_datafile');
        if oldversion  == 0
            data.atolls_datafile = 'none';
        end

%         oldversion = exist('invert.solve_for_fault_size');
%         if oldversion  == 0
%             invert.solve_for_fault_size = 'no';
%         end


%         oldversion = exist('invert.solve_for_InSAR_offset');
%         if oldversion  == 0
%             invert.solve_for_InSAR_offset = 'no';
%         end

        if strcmp(variable_patches_along_strike_with_depth, 'no') == 1
            visual_slip_mean = reshape(patch_mean, total_n_down_dip_patches, total_n_along_strike_patches);
        end

        % Find median
%         if strcmp(invert.solve_for_fault_size, 'yes') == 1
%             patch_median = zeros(n_slip_patches_on_each_fault_strand,1);
%             rake_median = zeros(n_slip_patches_on_each_fault_strand,1);
%             for i = 1:total_n_slip_patches
%                 patch_median(i,1) = median(slip_keep(i,onlyonkeptpatches(i,:)));
%                 rake_median(i,1) = median(rake_keep(i,onlyonkeptpatches(i,:)));
%             end
%             patch_median(isnan(patch_median)==1) = 0;
%             rake_median(isnan(rake_median)==1) = 0;
%         else
            patch_median = median(slip_keep, 2); % to take the median of each row
            rake_median = median(rake_keep,2);
        %end
        
        
        if strcmp(variable_patches_along_strike_with_depth, 'no') == 1
            visual_slip_median = reshape(patch_median, total_n_down_dip_patches, total_n_along_strike_patches);
        end

        % Find mode
         nbins = 30;
%         for i = 1: total_n_slip_patches
%             [s_counts(i,:), s_centers(i,:)] = hist(slip_keep(i,:), nbins);   % count the number in each bin and the center of each bin
%             [r_counts(i,:), r_centers(i,:)] = hist(rake_keep(i,:), nbins);
%         end
%             [M,s_I] = max(s_counts'); % If A is a matrix, then max(A) is a row vector containing the maximum value of each column
%             [M,r_I] = max(r_counts');
%         for i = 1: total_n_slip_patches
%             patch_mode(i,1) = s_centers(i,s_I(i));
%             rake_mode(i,1) = r_centers(i,r_I(i));
%         end

%         % Find mode 2d way
        for i = 1: total_n_slip_patches
%             if strcmp(invert.solve_for_fault_size, 'yes') == 1
%                slipkeeprow = slip_keep(i,onlyonkeptpatches(i,:));
%                rakekeeprow = rake_keep(i,onlyonkeptpatches(i,:));
%             else
                slipkeeprow = slip_keep(i,:);
                rakekeeprow = rake_keep(i,:);
            %end
            [N,Xedges,Yedges] = histcounts2(slipkeeprow,rakekeeprow,30);
            [~,I] = max(N(:));
            [I_row, I_col] = ind2sub(size(N),I);
            patch_mode(i,1) = Xedges(I_row);
            rake_mode(i,1) = Yedges(I_col);
        end
        
        if strcmp(variable_patches_along_strike_with_depth, 'no') == 1
            visual_slip_mode = reshape(patch_mode, total_n_down_dip_patches, total_n_along_strike_patches);
        end 

        % Find std
        if strcmp(invert.solve_for_fault_size, 'yes') == 1
            patch_std = zeros(n_slip_patches_on_each_fault_strand,1);
            rake_std = zeros(n_slip_patches_on_each_fault_strand,1);
            for i = 1:total_n_slip_patches
                patch_std(i,1) = std(slip_keep(i,onlyonkeptpatches(i,:)));
                rake_std(i,1) = std(rake_keep(i,onlyonkeptpatches(i,:)));
            end
        else
            rake_std = std(rake_keep.');
            patch_std = std(slip_keep.'); % calculates standard deviation of each column of slip_keep
        end
        if strcmp(variable_patches_along_strike_with_depth, 'no') == 1
            visual_std = reshape(patch_std, total_n_down_dip_patches, total_n_along_strike_patches);
        end
        %std_scaling = (1./visual_std)./(max(1./patch_std));
        std_scaling = ones(total_n_down_dip_patches, total_n_along_strike_patches);

        % Find maximum likelihood
        % %[~, I] = max(sum(logL_keep));      %test. this should equal highest_posterior: posterior_keep(I).    also sum since we need to sum the posterior for each patch
        % [~, I] = max(logL_keep);
        % patch_mostlikely = slip_keep(:, I);
        % % if strcmp(invert.smoothing, 'VK') == 1;
        % %     alpha_mostlikely = alpha_keep(1:n_fault_strands,I);
        % % end
        % rake_mostlikely = rake_keep(:,I);
        % visual_patch_mostlikely = reshape(patch_mostlikely, total_n_down_dip_patches, total_n_along_strike_patches);
        %WRMS_mostlikely = logL_keep(:,I);
        %WRMS_mostlikely = calc_WRMS(patch_mostlikely, G_mostlikely, d, inv_sigma_d);

        % Find maximum likelihood
        %if strcmp(invert.smoothing, 'none') ==1 && strcmp(invert.regularise_moment, 'no') ==1 && strcmp(invert.solve_for_dip, 'no') == 1; % posterior keep is only black if we JUST do a bayesian inversion with NOTHING ELSE
        oldversion = exist('logposterior_keep');
        if oldversion  == 0
            logposterior_keep = posterior_keep;
        end
        
        if strcmp(invert.smoothing, 'none') ==1 && strcmp(invert.regularise_moment, 'no') ==1 % posterior keep is only black if we JUST do a bayesian inversion with NOTHING ELSE
            [~, I] = max(logL_keep);
            patch_MAP = slip_keep(:, I);
        else
            [~, I] = max(sum(logposterior_keep,1));  % sum for each slip patch
            patch_MAP = slip_keep(:, I);   
        end
        
        oldversion = exist('alpha2_keep');
        if oldversion  == 0
            alpha2_keep = alpha_keep;
        end
        if strcmp(invert.smoothing, 'VK') == 1 || strcmp(invert.smoothing, 'laplacian') == 1
            alpha2_mostlikely = alpha2_keep(1:n_fault_strands_for_smoothing,I);
        end
        rake_mostlikely = rake_keep(:,I);
        if strcmp(variable_patches_along_strike_with_depth, 'no') == 1
            visual_patch_MAP = reshape(patch_MAP, total_n_down_dip_patches, total_n_along_strike_patches);
        end
        %WRMS_mostlikely = logL_keep(:,I);


        % Calculate best offset

        if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
           for i = 1:n_InSAR_scenes
                offset_mean(i,1) = mean(offset_keep(i,:));       % calculate the mean of the offset, but only the offsets relating to InSAR (the offsets relating to other datasets will be 0)
                [counts, centers] = hist(offset_keep(i,:), nbins);   % count the number in each bin and the center of each bin
                [M,I] = max(counts'); % If A is a matrix, then max(A) is a row vector containing the maximum value of each column
                offset_mode(i,1) = centers(1,I);
           end
           offset_MAP = offset_keep(:, I);
        else
            offset_mean = offset_initial;
            offset_MAP = offset_initial;
            offset_median = offset_initial;
            offset_mode = offset_initial;
        end
        
        % Calculate best ramp
        if strcmp(invert.solve_for_InSAR_ramp, 'yes') == 1
            for i = 1:(n_InSAR_scenes*3)
                ramp_mean(i,1) = mean(ramp_keep(i,:));       % calculate the mean of the ramp, but only the ramps relating to InSAR (the ramps relating to other datasets will be 0)
                [counts, centers] = hist(ramp_keep(i,:), nbins);   % count the number in each bin and the center of each bin
                [M,I] = max(counts'); % If A is a matrix, then max(A) is a row vector containing the maximum value of each column
                ramp_mode(i,1) = centers(1,I);
           end
           ramp_MAP = ramp_keep(:, I);
        else
            ramp_mean = ramp_initial;
            ramp_MAP = ramp_initial;
            ramp_median = ramp_initial;
            ramp_mode = ramp_initial;
        end

        % Calculate residuals to the different models

        if strcmp(invert.variable_rake, 'yes') == 1

            % mean    
                %theta = 180 - rake_mean;
                theta = degtorad(rake_mean);
                       for r = 1 : total_n_slip_patches             % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake.          
                                G_mean(:,r) = cos(theta(r)) * G_ss(:,r) + sin(theta(r)) * G_ds(:,r);    % I THINK. check yo' trigonometry    
                       end 
                %WRMS_mean = calc_WRMS(patch_mean, G_mean, d, inv_sigma_d, offset_mean);

            % median   
                %theta = 180 - rake_median;
                theta = degtorad(rake_median);        
                       for r = 1 : total_n_slip_patches             % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake.
                                G_median(:,r) = cos(theta(r)) * G_ss(:,r) + sin(theta(r)) * G_ds(:,r);    % I THINK. check yo' trigonometry 
                       end     
                %WRMS_median = calc_WRMS(patch_median, G_median, d, inv_sigma_d, offset_mmedian);

            % mode   
                %theta = 180 - rake_mode;
                theta = degtorad(rake_mode);    
                       for r = 1 : total_n_slip_patches             % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake.
                                 G_mode(:,r) = cos(theta(r)) * G_ss(:,r) + sin(theta(r)) * G_ds(:,r);    % I THINK. check yo' trigonometry 
                       end         
                %WRMS_mode = calc_WRMS(patch_mode, G_mode, d, inv_sigma_d, offset_mode);

            % most_likely   
                %theta = 180 - rake_mostlikely;
                theta = degtorad(rake_mostlikely);         
                       for r = 1 : total_n_slip_patches             % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake.
                                 G_mostlikely(:,r) = cos(theta(r)) * G_ss(:,r) + sin(theta(r)) * G_ds(:,r);    % I THINK. check yo' trigonometry  
                       end      
                %WRMS_mostlikely = calc_WRMS(patch_MAP, G_mostlikely, d, inv_sigma_d, offset_MAP);  

        elseif strcmp(invert.variable_rake, 'no') == 1
                %WRMS_mean = calc_WRMS(patch_mean, G, d, inv_sigma_d);
                %WRMS_median = calc_WRMS(patch_median, G, d, inv_sigma_d);
                %WRMS_mode = calc_WRMS(patch_mode, G, d, inv_sigma_d);
                %WRMS_mostlikely = calc_WRMS(patch_MAP, G, d, inv_sigma_d);
        end

        % Calculate 95% confidence intervals
        patch_conf_intervals = prctile(slip_keep',[2.5 97.5]);
        patch_conf_intervals = patch_conf_intervals'; %one line per slip patch

        ninetyfive_percent_conf = diff(patch_conf_intervals,1,2);

end

%% Plot confidence intervals

%can only use this if you used a version of slipBERi that calculated CI_low and CI_high
% plot iterations vs confidence intervals and change in confidence intervals
%     figure
%     plot( store_number_at_sens_test, CI_low', 'b');
%     hold on;
%     plot( store_number_at_sens_test, CI_high', 'r');
%     title('95% confidence');
%     legend('Low', 'High');
%     xlabel('Iterations')
%     ylabel('Slip (meters)');
% 
%     figure
%     plot( store_number_at_sens_test, diff(CI_low)', 'b');
%     hold on;
%     plot( store_number_at_sens_test, diff(CI_high)', 'r');
%     title('diff in 95% confidence');
%     legend('Low', 'High');
%     xlabel('Iterations')
%     ylabel('Slip (meters)');
%     
%     figure
%     plot( store_number_at_sens_test, mean(diff(CI_low))', 'b');
%     hold on;
%     plot( store_number_at_sens_test, mean(diff(CI_high))', 'r');
%     title('mean diff in 95% confidence');
%     legend('Low', 'High');
%     xlabel('Iterations')
%     ylabel('Slip (meters)');


% % restrospectively
% 
% sens_test_for_conf = sens_test;
% sens_test_for_conf(sens_test_for_conf > size(slip_keep,2)) =[];
% sens_test_for_conf(end) = []; %coz that's nan
% 
% for i = 1: length(sens_test_for_conf)
% 
%      slip_to_calc_conf = slip_keep(:, 1:sens_test_for_conf(i));
% %     n_entries = size(slip_to_calc_conf,2);
% %     SEM_beginning_to_end(:,i) = (std( slip_to_calc_conf') ./sqrt(n_entries))';               % Standard Error
% %     ts = tinv([0.05  0.95],(n_entries-1));      % T-Score
% %     CI_low_beginning_to_end(:,i) = mean(slip_to_calc_conf')' + (ts(1)*SEM_beginning_to_end(:,i));                     % Confidence Intervals
% %     CI_high_beginning_to_end(:,i) = mean(slip_to_calc_conf')' + ts(2)*SEM_beginning_to_end(:,i);                      % Confidence Intervals
% %     clear slip_to_calc_conf
%     confidence_intervals_beginning_to_end(:,:,i) = prctile(slip_to_calc_conf',[5 95]);
%     
% end
% 
% CI_low_beginning_to_end(:,:) = confidence_intervals_beginning_to_end(1,:,:);
% CI_high_beginning_to_end(:,:) = confidence_intervals_beginning_to_end(2,:,:);

%confidence_intervals_diff = diff(confidence_intervals_beginning_to_end);

% % plot confidence interval for all slip patches
% figure('position', [100, 300, 800, 1000])
% subplot(3,1,1)
% plot( sens_test_for_conf, CI_low_beginning_to_end', 'b');
% hold on;
% plot( sens_test_for_conf, CI_high_beginning_to_end', 'r');
% title('95% confidence');
% ylim([0 max(CI_low_beginning_to_end(:,1))])
% %legend('Low', 'High');
% %xlabel('Iterations')
% %ylabel('Slip (meters)');
% 
% % plot difference in confidence interval for all slip patches
% %figure
% subplot(3,1,2)
% plot( sens_test_for_conf, diff(CI_low_beginning_to_end)', 'b');
% hold on;
% plot( sens_test_for_conf, diff(CI_high_beginning_to_end)', 'r');
% title('95% confidence difference');
% ylim([min(diff(CI_low_beginning_to_end(:,1)')) max(diff(CI_low_beginning_to_end(:,1)'))])
% %legend('Low', 'High');
% %xlabel('Iterations')
% ylabel('Slip (meters)');
% 
% % this is the mean of all of the patches (coz I was getting may be 10 outliers, but for all I know they could be poorly constrained ones at hte bottom.
% %figure
% subplot(3,1,3)
% plot( sens_test_for_conf, mean(diff(CI_low_beginning_to_end))', 'b');
% hold on;
% plot( sens_test_for_conf, mean(diff(CI_high_beginning_to_end))', 'r');
% title('95% confidence mean difference');
% %legend('Low', 'High');
% xlabel('Iterations')
% %ylabel('Slip (meters)');

%% If in testing mode, plot 95% conf with true value

% see whether the true values in testing mode plot with in the 95%
if strcmp(testing.testing_mode, 'yes') && strcmp(invert.smoothing, 'tikhonov') == 0
   
    plot_numbers = 1: total_n_slip_patches;
    
%     % if you want to plot patches 1:end in ROWS instead of columns
%     % (default) then uncomment this
     A = 1: total_n_slip_patches;
     A = reshape(A, total_n_down_dip_patches, n_along_strike_patches);
     A = fliplr(A);
     A = A';
     plot_numbers = reshape(A, total_n_slip_patches, 1);
%      A = 1: total_n_slip_patches;
%      A = reshape(A, total_n_down_dip_patches, total_n_along_strike_patches);
%       A = A';    
%       B=flipud(A);
%     %A = fliplr(A);
%      A = A';
%     A = flipud(A);
%     A = A';
%     plot_numbers = reshape(B, total_n_slip_patches, 1);
%     %plot_numbers2 = reshape(A, total_n_slip_patches, 1);

     %plot_numbers2 = reshape(flipud(A), total_n_slip_patches, 1);
    
    plot_numbers2 = plot_numbers;

    figure('position', [100, 300, 1800, 1000])
    hold all
    %plot( 1:100, patch_conf_intervals(plot_numbers,:), 'LineWidth', 1.5);
    X = [1:total_n_slip_patches; 1:total_n_slip_patches];
    Y = [patch_conf_intervals(plot_numbers,1), patch_conf_intervals(plot_numbers,2)]';
    plot( X, Y, 'r', 'LineWidth', 1.5);
    hold on;
    %scatter( 1:total_n_slip_patches, patch_mode(plot_numbers), 30, 'k', 'filled');
    scatter( 1:total_n_slip_patches, patch_mean(plot_numbers), 30, 'k', 'filled');
    xlabel('Patch number (Along strike)')
    ylabel('Slip (m)')
    title('95% confidence')
     
    if strcmp(testing.testing_mode, 'yes') == 1
       
        % true value
        load(testing.making_model, 'synthetic_slip');
        
        % how to find which patch number it is if the model and solution have a different number of slip patches?
        

%         if length(synthetic_slip) ~= length(patch_mode);
%             disp('I don''t have a better way of doing this yet');
%             prompt = 'Type the number slip patch that each true slip patch equates to on your own fault as a matrix';
%         truepatchplot_numbers = input(prompt);
%         else
            truepatchplot_numbers = 1 :length(synthetic_slip);
%         end
%         
%         % plot on top in blue
         scatter( truepatchplot_numbers, synthetic_slip(plot_numbers2), 200, 'gx', 'LineWidth', 3);
        
    end
     
%     legend('95% conf high', '95% conf low', 'Mode', 'True value')
%     legend('95% conf high', 'Mode', 'True value')
    
    % add dotty lines to separate rows
    %betweenrows =(total_n_along_strike_patches:total_n_along_strike_patches:total_n_slip_patches)+0.5;      % this was back when n_along_strike_patches was one value, intead of one value per row
    betweenrows = total_n_along_strike_patches +0.5;
    %Y = [0 ceil(max(max(patch_conf_intervals)))];
    Y = [0 6];
    Y = repmat(Y', 1, (total_n_down_dip_patches(1)));
    %X = [betweenrows ; betweenrows];
    %plot(X,Y, '--b');
    
    for i = 1 : total_n_down_dip_patches-1               % do the loop to plot all the dotty lines between allh te rows
        X = [betweenrows + (total_n_along_strike_patches*(i-1)); betweenrows+ (total_n_along_strike_patches*(i-1))];
        plot(X,Y, '--b');
    end
    
    
    xlim([0 total_n_slip_patches])
    
    for i = 1 : total_n_down_dip_patches
        x1 = (total_n_along_strike_patches/2)+(total_n_along_strike_patches*(i-1))-2;
        y1 = 2.5;
        txt1 = ['Row ', num2str(i)];
        %text(x1,y1,txt1, 'Fontsize', 22)
    end
    
    if strcmp(testing.testing_mode, 'yes') == 1
        legend('95% conf', 'Patch mode', 'True value')
    else
        legend('95% conf', 'Patch mode')
    end
    
    ylim([0 6])
    
end

%% Plot the result, mean, flat

if strcmp(display.plotmean, 'yes') + strcmp(display.plotmode, 'yes') + strcmp(display.plotmedian, 'yes') ~=0

scaling_factor = 0.5;

if strcmp(variable_patches_along_strike_with_depth, 'no') == 1;         % can't plot using imagesc if there are a different number of patches along strike

    
            order = flipud(spatial_model1);
            order = reshape(order, 1, []);
            order = order';
    
    
        things_to_get_around_to_plotting = {};
        visual_things_to_get_around_to_plotting = {};
        rake_to_get_around_to_plotting ={};
        if strcmp(display.plotmean, 'yes') == 1;
             things_to_get_around_to_plotting = {'patch_mean'};
             visual_things_to_get_around_to_plotting = {'visual_slip_mean'};
             rake_to_get_around_to_plotting = {'rake_mean'};
        end
        if strcmp(display.plotmode, 'yes') == 1;
            things_to_get_around_to_plotting{end+1,1} = {'patch_mode'};
            visual_things_to_get_around_to_plotting{end+1,1}  = {'visual_slip_mode'};
            rake_to_get_around_to_plotting{end+1,1} = {'rake_mode'};
        end
        if strcmp(display.plotmedian, 'yes') == 1;
            things_to_get_around_to_plotting{end+1,1} = {'patch_median'};
            visual_things_to_get_around_to_plotting{end+1,1}  = {'visual_slip_median'};
            rake_to_get_around_to_plotting{end+1,1} = { 'rake_median'};
        end
        if strcmp(display.plotMAP, 'yes') == 1;
            things_to_get_around_to_plotting{end+1,1} = {'patch_MAP'};
            visual_things_to_get_around_to_plotting{end+1,1}  = {'visual_patch_MAP'};
            rake_to_get_around_to_plotting{end+1,1} = {'rake_mostlikely'};
        end        

    
        for i = 1: size(visual_things_to_get_around_to_plotting,1)     % loop through all the things we're plotting
        
                thing_plotting_atm = things_to_get_around_to_plotting{i};
                visual_thing_plotting_atm = visual_things_to_get_around_to_plotting{i};
                rake_plotting_atm = rake_to_get_around_to_plotting{i};

                if iscell(visual_thing_plotting_atm)
                    visual_thing_plotting_atm = cell2mat(visual_thing_plotting_atm);
                end
                if iscell(thing_plotting_atm)
                    thing_plotting_atm = cell2mat(thing_plotting_atm);
                end
                if iscell(rake_plotting_atm)
                    rake_plotting_atm = cell2mat(rake_plotting_atm);
                end
                

                %std_scaling = (1./visual_std)./(max(1./patch_std));     
                %std_scaling(visual_std>visual_slip_mean) = 0.01;          % set any slip patches to zero if std of a patch is greater than the slip value
                std_scaling = ones(total_n_down_dip_patches, total_n_along_strike_patches);

                % plot mode and std
                figure('position', [100, 300, 600, 800])
                subplot(2,1,1)
                h = eval(visual_thing_plotting_atm);
                h = reshape( flipud(h(order)), total_n_down_dip_patches,[]);
                %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_mode);
                imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2],h, 'AlphaData', std_scaling);
                %imagesc([spatial_model3(1)/1000/2, fault_length(1)/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_mode);
                axis equal tight
                hot2 = hot;
                hot2(1:8,:) = [];   % making the top colour red, not v dark red
                hot2 = [hot2; ones(5,3)];  % adding more white to the bottom
                colorbar
                colormap(hot2)
                colormap(flipud(colormap))
                %xlabel('Along Strike'); 
                ylabel('Down dip'); ylabel(colorbar, 'Slip (m)');
                title(thing_plotting_atm)
                %str = ['WRMS = ', num2str(WRMS_mode)];
                %text(-2, -2, str);

                
                    if strcmp(invert.variable_rake, 'yes') == 1;
                        [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( eval(rake_plotting_atm), eval(thing_plotting_atm), spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model );
                    elseif strcmp(invert.variable_rake, 'no') == 1;
                        [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( disloc_model(5,:), eval(thing_plotting_atm), spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model );   
                    end
                    subplot(2,1,1); hold on;
                    %quiver(x,y,u_arrow,v_arrow, scaling_factor, 'k');

                ax2 = subplot(2,1,2);
                imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2],reshape( flipud(visual_std(order)), total_n_down_dip_patches,[]))
                axis equal tight
                colorbar
                colormap( ax2,'jet')
                xlabel('Along Strike (km)'); %ylabel('Down dip'); 
                ylabel(colorbar, 'std (m)');
                title('Standard deviation')
                
                clear thing_plotting_atm;
                clear visual_thing_plotting_atm;
                clear rake_plotting_atm;
                
            
        end
                

end

end

%% Compare all slips


if strcmp(display.plotallslips, 'yes') == 1

    %colormap('hot')
    
    scaling_factor = 0.35;
    
            order = flipud(spatial_model1);
            %order = fliplr(order);
            order = reshape(order, 1, []);
            order = order';
            
            %order2=order;
            order2=order;
            order2 = flipud(order2);
            order2 = reshape(order2, 1, []);
           
    
    % first compare slip distributions ****************************************
        figure('position', [300, 300, 1600, 1000])

        if strcmp(testing.testing_mode, 'yes') == 1
            load(testing.making_model, 'synthetic_slip');
            if size(synthetic_slip,2) > size(synthetic_slip,1) ==1
                synthetic_slip = synthetic_slip';
            end
            cmaxvalue = max([patch_mean; patch_mode; patch_median; patch_MAP; synthetic_slip]);
        else
            cmaxvalue = max([patch_mean; patch_mode; patch_median; patch_MAP]);
            %cmaxvalue = max([patch_mean; patch_mode; patch_median]);
        end
        
        if strcmp(testing.testing_mode, 'yes') == 1
            % actual slip distribution, only use if testing
            load(testing.making_model, 'visual_slip', 'synthetic_slip')
            subplot(3,1,1);
            %imagesc([0, sum(fault_length)/1000],[0, (fault_width(1))/1000],visual_slip); hold on;
            %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2],visual_slip);
            imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2],visual_slip);
            %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip); hold on;
            hold on;
            axis equal tight; colorbar;
            ylim([0, (fault_width(1))/1000])
            xlim([0, sum(fault_length)/1000])
            
            true_model = load(testing.making_model, 'synthetic_disloc_model');
            rake_true = true_model.synthetic_disloc_model(5,:); 
            spatial_model2_for_true = spatial_model2;
            n_down_dip_patches_for_plotting = total_n_down_dip_patches;
            n_along_strike_patches_for_plotting = total_n_along_strike_patches;  
            n_down_dip_patches_for_true = total_n_down_dip_patches(1);
            n_along_strike_patches_for_true = sum(n_along_strike_patches);
            %[ x,y,u_arrow,v_arrow ] = prepare_plot_rake( rake_true, synthetic_slip, spatial_model2_for_true, n_down_dip_patches_for_true, n_along_strike_patches_for_true, true_model.synthetic_disloc_model );
            true5 =load(testing.making_model, 'total_n_along_strike_patches'); 
            true6 =load(testing.making_model, 'total_n_down_dip_patches'); 
            true7 =load(testing.making_model, 'spatial_model3'); 
            [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( rake_true, synthetic_slip, true1.spatial_model2, true7.spatial_model3, true6.total_n_down_dip_patches, true5.total_n_along_strike_patches, true_model.synthetic_disloc_model );
            quiver(x,y,u_arrow,v_arrow, scaling_factor, 'k');            
    
            xlabel('along strike')
            ylabel('down dip')
            title('Actual solution')          
            %caxis([0, max([patch_mean; patch_mode; patch_median; patch_MAP])]);
            caxis([0, cmaxvalue]);
                    %caxis([0, 25]);
        end 
        


        % mean
        if strcmp(testing.testing_mode, 'yes') == 1
            subplot(3,2,3)
        else
            subplot(2,2,1)
        end
        %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_mean); hold on;
        imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2], reshape( flipud(patch_mean(order)), total_n_down_dip_patches,[]) );
        hold on; 
        %imagesc([spatial_model3(1)/1000/2, fault_length(1)/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_mean);hold on;
        axis equal tight
        colorbar          
            if strcmp(invert.variable_rake, 'yes') == 1;
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( rake_mean, patch_mean, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model );
            elseif strcmp(invert.variable_rake, 'no') == 1;
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( disloc_model(5,:), patch_mean, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model );   
            end
            quiver(x(order2),y(order2),u_arrow,v_arrow, scaling_factor, 'k');
                     
        %xlabel('Along Strike (km)'); 
        ylabel('Down dip (km)'); %ylabel(colorbar, 'Slip (m)');
        title('Imaginary EQ slip, with VK constraint')
        title('Mean')
        %caxis([0, max([patch_mean; patch_mode; patch_median; patch_MAP])]);   
        %caxis([0, max([patch_mean; patch_mode; patch_median])]);
        caxis([0, cmaxvalue]);
                %caxis([0, 25]);
        %str = ['WRMS = ', num2str(WRMS_mean)];
        %text(-2, -2, str);

        % mode
        if strcmp(testing.testing_mode, 'yes') == 1
           subplot(3,2,4)
        else
            subplot(2,2,2)
        end
        %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_mode); hold on;
        imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2], reshape( flipud(patch_mode(order)), total_n_down_dip_patches,[]));
        hold on; 
        %imagesc([spatial_model3(1)/1000/2, fault_length(1)/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_mode);hold on;
        axis equal tight
        colorbar
            if strcmp(invert.variable_rake, 'yes') == 1;
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( rake_mode, patch_mode, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model );
            elseif strcmp(invert.variable_rake, 'no') == 1;
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( disloc_model(5,:), patch_mode, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model );   
            end
            quiver(x(order2),y(order2),u_arrow,v_arrow, scaling_factor, 'k');
        %xlabel('Along Strike'); ylabel('Down dip'); 
        ylabel(colorbar, 'Slip (m)');
        title('Imaginary EQ slip, with VK constraint')
        title('Mode')
        %caxis([0, max([patch_mean; patch_mode; patch_median; patch_MAP])]); 
        %caxis([0, max([patch_mean; patch_mode; patch_median])]);
        caxis([0, cmaxvalue]);
                %caxis([0, 25]);
        %str = ['WRMS = ', num2str(WRMS_mode)];
        %text(-2, -2, str);

        % median
        if strcmp(testing.testing_mode, 'yes') == 1
           subplot(3,2,5)
        else
            subplot(2,2,3)
        end
        %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_median); hold on;
        imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2], reshape( flipud(patch_median(order)), total_n_down_dip_patches,[]));
        hold on; 
        %imagesc([spatial_model3(1)/1000/2, fault_length(1)/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_slip_median);hold on;       
        axis equal tight
        colorbar
            if strcmp(invert.variable_rake, 'yes') == 1
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( rake_median, patch_median, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches,disloc_model );
            elseif strcmp(invert.variable_rake, 'no') == 1
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( disloc_model(5,:), patch_median, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model );   
            end
            quiver(x(order2),y(order2),u_arrow,v_arrow, scaling_factor, 'k');
        %xlabel('Along Strike'); ylabel('Down dip'); 
        ylabel(colorbar, 'Slip (m)');
        title('Imaginary EQ slip, with VK constraint')
        title('Median')
        %caxis([0, max([patch_mean; patch_mode; patch_median; patch_MAP])]);
        %caxis([0, max([patch_mean; patch_mode; patch_median])])
        caxis([0, cmaxvalue]);
                %caxis([0, 25]);
        %str = ['WRMS = ', num2str(WRMS_median)];
        %text(-2, -2, str);
        
        % MAP
        if strcmp(testing.testing_mode, 'yes') == 1
           subplot(3,2,6)
        else
            subplot(2,2,4)
        end
        %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_patch_MAP); hold on;
        imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2], reshape( flipud(patch_MAP(order)), total_n_down_dip_patches,[]));
        %imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2],visual_patch_MAP);
        hold on; 
        %imagesc([spatial_model3(1)/1000/2, fault_length(1)/1000], [spatial_model3(1)/1000/2, (fault_width(1)-(spatial_model2(1)/2))/1000],visual_patch_mostlikely);       
        axis equal tight
        colorbar
            if strcmp(invert.variable_rake, 'yes') == 1;
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( rake_mostlikely, patch_MAP, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, disloc_model);
            elseif strcmp(invert.variable_rake, 'no') == 1;
                [ x,y,u_arrow,v_arrow ] = prepare_plot_rake( disloc_model(5,:), patch_MAP, spatial_model2, spatial_model3, total_n_down_dip_patches, total_n_along_strike_patches, testing, disloc_model );   
            end
            quiver(x(order2),y(order2),u_arrow,v_arrow, scaling_factor, 'k');
        xlabel('Along Strike'); %ylabel('Down dip'); 
        ylabel(colorbar, 'Slip (m)');
        title('Imaginary EQ slip, with VK constraint')
        title('MAP')
        %caxis([0, max([patch_mean; patch_mode; patch_median; patch_MAP])]); 
        %caxis([0, max([patch_mean; patch_mode; patch_median])]);
        caxis([0, cmaxvalue]);
        %caxis([0, 25]); 
        %str = ['WRMS = ', num2str(WRMS_mostlikely)];
        %text(-2, -2, str);
        
        %colormap('hot')
        %colormap(flipud(colormap))
%         hot2 = hot;
%         hot2(1:8,:) = [];      % chopping off start of colorbar to get rid of black - i want the maximum colour to be dark red
%         colormap(hot2)
%         colormap(flipud(colormap))
        hot2 = hot;
        hot2(1:5,:) = [];      % chopping off start of colorbar to get rid of black - i want the maximum colour to be dark red  
        hot2 = [hot2; ones(5,3)];  % adding more white to the bottom
        colormap(hot2)
        colormap(flipud(colormap))
        
end



%% Probability with time

if strcmp(display.plotprob, 'yes')
    % Have a look and see if you escaped the burn-in period
    figure
    %posterior_total = sum(posterior_keep,1);                  % this sums the probability, if there are multiple strands, and if there's just one strand it does nothing
    plot( 1:length(logposterior_keep), logposterior_keep)
    %axis ij % If you'd like to flip the axis - maybe do this if you're taking the log of the probabilities
    xlabel('Iterations')
    ylabel('Unscaled Probability')
end

%%  Patch on-off-edness

if strcmp(invert.solve_for_fault_size, 'yes') ==1

    onoffidentifyerkeeppercentage = sum(onlyonkeptpatches,2) / (store_number-burn_in_remove_number);
    visualonoffidentifyerkeep = reshape(onoffidentifyerkeeppercentage, total_n_down_dip_patches, n_along_strike_patches);
    figure('position', [600, 300, 800, 600])
    imagesc([spatial_model3(1)/1000/2, (sum(fault_length)-(spatial_model2(1)/2))/1000], [spatial_model3(1)/1000/2, (fault_width(1)/1000/sin(deg2rad(disloc_model(4,1))))-spatial_model3(1)/1000/2],visualonoffidentifyerkeep)
    title('Patch on-off-edness (%)');
    h = colorbar;
    colormap(hot)
    colormap(flipud(colormap))
    caxis([0 1])

    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    
%     h = colorbar;
%     set( h, 'YDir', 'reverse' );
    
    % work out how many times a patch has to be on to be significant
    %binoprob = binopdf(onoffidentifyerstoretotal,store_number,0.2);     % Y = binopdf(X,N,P);

    
    
end


%% Plot histograms for some slip patches

if strcmp(display.plothists, 'plothistsample') == 1
        % Select some slip patches randomly
        a = 1;
        b = total_n_slip_patches;
        random = (b-a).*rand(5,1) + a;
        random = round(random);
        random = sort(random);

        % Show which patches I'm actually looking at
        figure('position', [100, 350, 1600, 1200])
        subplot(3,2,1)
        location_matrix = zeros( total_n_down_dip_patches, total_n_along_strike_patches);
        location_matrix(random(1)) = 1;
        location_matrix(random(2)) = 1;
        location_matrix(random(3)) = 1;
        location_matrix(random(4)) = 1;
        location_matrix(random(5)) = 1;
        imagesc(location_matrix)
        xlabel('Along strike')
        ylabel('Down dip')
        %ax1 = subplot(3,2,1);
        %colours = [ 1 1 1; 1 0 0];
        %colormap(ax1, colours);
        %caxis([0 1])
        title('Location of randomly sampled slip patches')

        % Plot histograms for each of these patches
        for i = 2:6

        subplot(3,2,i)
        colormap default
        hist( slip_keep( random(i-1), :),nbins)
        hold on
        y = [0 max(M)];
        x_mean = [patch_mean(random(i-1)) patch_mean(random(i-1))];
        x_mode = [patch_mode(random(i-1)) patch_mode(random(i-1))];
        x_median = [patch_median(i) patch_median(i)];
        plot(x_mean,y,'r');
        plot(x_mode,y,'g');
        plot(x_median,y,'k');
        plot(patch_conf_intervals(i,1));
        plot(patch_conf_intervals(i,2));
           if strcmp(testing.testing_mode, 'yes')
              load(testing.making_model, 'synthetic_slip');
              x_true = synthetic_slip(i); 
              plot(x_true, y, 'm');
           end
        xlabel('Slip (m)')
        ylabel('Number of occurences')
        title(['Slip patch ', num2str(random(i-1))])
        ylim([0, iterations/3]);

        end

        % Put a legend, but only on the second histogram
        subplot(3,2,2)
        legend('Sampling', 'Mean', 'Mode');
        legend('location', 'northeast')
            if strcmp(testing.testing_mode, 'yes')
                legend('Sampling', 'Mean', 'Mode', 'Median', 'True');
            else
                legend('Sampling', 'Mean', 'Mode', 'Median','95% Confidence intervals');
            end
    legend('location', 'northeast')
        
        
end



%% Surface displacment arrows
% can only do this for testing, because for InSAR we only have LOS

    if strcmp(display.plotsurfacedisp, 'yes') + strcmp(testing.testing_mode, 'yes') == 2
       if strcmp(testing.testing_mode, 'yes') ==1
            load(testing.making_model, 'u');
            from_making_model = load(testing.making_model, 'locs');
           
            scale = 340; 
            
            figure('position', [100, 300, 800, 800])
            ax1 = subplot(2,1,1);
            
            % plot fault trace       
            for j = 1 : n_fault_strands
                x = [fault_coords(j,1), fault_coords(j,3)];
                y = [fault_coords(j,2), fault_coords(j,4)];
                hold on;
                plot(x, y, 'm', 'Linewidth', 2)
            end

            % true data
            displ_actual_E = u(1,:);
            displ_actual_N = u(2,:);
            quiver( from_making_model.locs(1,:), from_making_model.locs(2,:), displ_actual_E* scale, displ_actual_N* scale, 'b', 'Autoscale', 'off');
            hold on;


            % model prediction       
            displ_pred_E = G_E*patch_MAP;         % d = Gm
            displ_pred_N = G_N*patch_MAP;
            quiver( locs(1,1: length(displ_actual_E)), locs(2,1: length(displ_actual_E)), displ_pred_E' * scale, displ_pred_N' * scale, 'g', 'Autoscale', 'off');
            legend('Fault', 'Actual data', 'Model prediction', 'Location', 'northeastoutside');
            load(testing.making_model, 'faultx', 'faulty');
            plot( faultx/1000, faulty/1000, 'k');
            xlabel('x'); ylabel('y');
            %daspect([1 1 1])
            title('Maximum likelihood solution');
            axis tight




            % let's look at residuals, too ********************************************
            ax2 = subplot(2,1,2);
            
            resid_E = displ_actual_E' - displ_pred_E;
            resid_N = displ_actual_N' - displ_pred_N;
            for j = 1 : n_fault_strands
                x = [fault_coords(j,1), fault_coords(j,3)];
                y = [fault_coords(j,2), fault_coords(j,4)];
                hold on;
                plot(x, y, 'm', 'Linewidth', 2)
            end
            quiver( locs(1,1:length(resid_E)), locs(2,1:length(resid_N)), resid_E'* scale, resid_N'* scale, 'r', 'Autoscale', 'off')
            legend('Fault', 'Residuals', 'Location', 'northeastoutside');
            xlabel('x'); ylabel('y');
            hold on;
            plot( faultx, faulty, 'k');
            %daspect([1 1 1])
       end
        
       linkaxes([ax1, ax2], 'xy')
       
    end
    
%% Compare data (InSAR/GPS) and solution

if strcmp(display.plotsurfacedisp, 'yes') + strcmp(testing.testing_mode, 'no') == 2;
   
    figure('position', [100, 350, 1800, 600])
    counter=1;
    
    for p= 1: max([n_InSAR_scenes 1])
    
        ax1 = subplot(max([n_InSAR_scenes 1]),3,counter);

        % First, plot InSAR, as long as there is some
        if strcmp(data.InSAR_datafile, 'none') ~= 1
            scatter( locs_InSAR(1,InSAR_identifyer==p), locs_InSAR(2,InSAR_identifyer==p), [], d_InSAR(InSAR_identifyer==p), 'filled')
            axis equal
            colorbar
        end

        scale_factor = 50000;
        % Then plot GPS, as long as there is some
        if strcmp(data.GPS_datafile_2d, 'none') == 0 || strcmp(data.GPS_datafile_3d, 'none') == 0
               hold on;
               quiver( locs_GPS_unique(1,:), locs_GPS_unique(2,:), d_GPS_e*scale_factor, d_GPS_n*scale_factor, 'k', 'Autoscale', 'off');   % i.e. plot every third row of the locations against E and N...
               ylabel('UTM y')
        end

        if strcmp(data.atolls_datafile, 'none') ~= 1 
           scatter( locs_atolls(1,:), locs_atolls(2,:), [], d_atolls, 'filled')
           hold on;
        end

        % plot fault trace       
            for j = 1 : n_fault_strands
                x = [fault_coords(j,1), fault_coords(j,3)];
                y = [fault_coords(j,2), fault_coords(j,4)];
                hold on;
                plot(x, y, 'm', 'Linewidth', 2)
            end
        
        if counter ==1
             title('Data');
        end
        
        counter = counter+1;


        % Second, plot surface displacement predicted from slip solution
        ax2 = subplot(max([n_InSAR_scenes 1]),3,counter);      
        
        if strcmp(invert.variable_rake, 'yes') == 1
            G = G_mostlikely;
         else
            G = G;
        end  

        d_hat =  G * patch_MAP;
%             model_disloc_model = disloc_model;
%             model_disloc_model(5,:) = rake_mostlikely;
%             model_disloc_model(6,:) = patch_MAP;
%             [u, flag]=disloc3d3(model_disloc_model, locs_InSAR, elastic_params.lambda, elastic_params.mu_okada);       % for each iteration, this gives sums displacement at each data location due to all slip patches
%             d_hat = sum(u.*los_vector_InSAR);
%         d_hat = d_hat';
%         d_hat = [d_hat; zeros(length(d_GPS),1)];

 
         if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
             d_hat(InSAR_identifyer==p) = d_hat(InSAR_identifyer==p) + offset_MAP(p);
         elseif strcmp(invert.solve_for_InSAR_ramp, 'yes') == 1
             rampvalues = Gramp*ramp_mode;
             rampvalues((end+1):length(d))=0; % GPS
             d_hat(InSAR_identifyer==p) = d_hat(InSAR_identifyer==p) + rampvalues(InSAR_identifyer==p,1);
         elseif strcmp(invert.smoothing, 'tikhonov') ==1
             d_hat(d_identifyer==1) = d_hat(d_identifyer==1) + best_offset_tik;
         end

        hold on;

        if strcmp(data.InSAR_datafile, 'none') ~= 1    
            scatter( locs_InSAR(1,InSAR_identifyer==p), locs_InSAR(2,InSAR_identifyer==p), [], d_hat(InSAR_identifyer==p), 'filled')
            hold on;
        end

        if strcmp(data.GPS_datafile_2d, 'none') + strcmp(data.GPS_datafile_3d, 'none') == 0
            d_hat_e_2d = d_hat((length(locs_InSAR)+1):2:(length(locs_InSAR) + length(d_GPS_2d)));
            d_hat_n_2d = d_hat((length(locs_InSAR)+2):2:(length(locs_InSAR) +length(d_GPS_2d)));
            d_hat_e_3d = d_hat(((length(locs_InSAR) + length(locs_GPS_2d)+1)):3:end);
            d_hat_n_3d = d_hat(((length(locs_InSAR) + length(locs_GPS_2d)+2)):3:end);
            d_hat_e = [d_hat_e_2d; d_hat_e_3d];
            d_hat_n = [d_hat_n_2d; d_hat_n_3d];
        elseif strcmp(data.GPS_datafile_3d, 'none') == 0   % just 3d 
            d_hat_e = d_hat((length(locs_InSAR)+1):3:end);
            d_hat_n = d_hat((length(locs_InSAR)+2):3:end);
        elseif strcmp(data.GPS_datafile_2d, 'none') == 0   % just 2d 
            d_hat_e = d_hat((length(locs_InSAR)+1):2:end);
            d_hat_n = d_hat((length(locs_InSAR)+2):2:end);
        end

        if strcmp(data.GPS_datafile, 'yes') == 1
            quiver( locs_GPS_unique(1,:), locs_GPS_unique(2,:), d_hat_e'*scale_factor, d_hat_n'*scale_factor, 'k', 'Autoscale', 'off');   % 2D so plot every second row of the locations against E and N...
        end


        if strcmp(data.atolls_datafile, 'none') ~= 1 
           scatter( locs_atolls(1,:), locs_atolls(2,:), [], d_hat((length(d_InSAR)+length(d_GPS)+1):end), 'filled')
           hold on;
        end


        axis equal
        colorbar

        % plot fault trace
            for j = 1 : n_fault_strands 
                x = [fault_coords(j,1), fault_coords(j,3)];
                y = [fault_coords(j,2), fault_coords(j,4)];
                hold on;
                plot(x, y, 'm', 'Linewidth', 2)

            end 

        if counter ==2
             title('Model (MAP)');
        end
        counter = counter+1;

        % Third, plot residuals
        ax3 = subplot(max([n_InSAR_scenes 1]),3,counter);
        residuals = (d - d_hat);%.^2;        % the first half are InSAR residuals, if we have InSAR. The second half are GPS residuals, if we have GPS.

        if strcmp(data.InSAR_datafile, 'none') ~= 1     % if we have InSAR data, then the first (length(locs_InSAR)) entries of d are for InSAR
            scatter( locs_InSAR(1,InSAR_identifyer==p), locs_InSAR(2,InSAR_identifyer==p), [], residuals(InSAR_identifyer==p), 'filled');
            hold on;
        end

        if strcmp(data.atolls_datafile, 'none') ~= 1     % if we have InSAR data, then the first (length(locs_InSAR)) entries of d are for InSAR
            scatter( locs_atolls(1,:), locs_atolls(2,:), [], residuals((length(d_InSAR)+length(d_GPS)+1):end), 'filled');
            hold on;
        end


        if strcmp(data.GPS_datafile_2d, 'none') + strcmp(data.GPS_datafile_3d, 'none') == 0
               %subplot(1,3,3); 
               hold on;
               % 2d
               quiver( locs_GPS(1, 1:2: length(locs_GPS_2d)), locs_GPS(2, 1:2:length(locs_GPS_2d)), residuals((length(locs_InSAR)+1):2:(length(locs_InSAR) + length(d_GPS_2d)))'*scale_factor, residuals((length(locs_InSAR)+2):2:(length(locs_InSAR) +length(d_GPS_2d)))'*scale_factor, 'k', 'Autoscale', 'off');   % 2D so plot every second row of the locations against E and N...
               % then 3d
               quiver( locs_GPS(1, (length(locs_GPS_2d)+1):3:end), locs_GPS(2, (length(locs_GPS_2d)+2):3:end), residuals(((length(locs_InSAR) + length(locs_GPS_2d)+1)):3:end)'*scale_factor, residuals(((length(locs_InSAR) + length(locs_GPS_2d)+2)):3:end)'*scale_factor, 'k', 'Autoscale', 'off');   % i.e. 3D so plot every third row of the locations against E and N...        
        elseif strcmp(data.GPS_datafile_2d, 'none') == 0; 
               %subplot(1,3,3); 
               hold on;
               quiver( locs_GPS(1,1:2: length(d_GPS)), locs_GPS(2,1:2: length(d_GPS)), residuals((length(locs_InSAR)+1):2:end)'*scale_factor, residuals((length(locs_InSAR)+2):2:end)'*scale_factor, 'k', 'Autoscale', 'off');   % i.e. plot every third row of the locations against E and N...              
        elseif strcmp(data.GPS_datafile_3d, 'none') == 0;
               %subplot(1,3,3); 
               hold on;
               quiver( locs_GPS(1,1:3: length(d_GPS)), locs_GPS(2,1:3: length(d_GPS)), residuals((length(locs_InSAR)+1):3:end)'*scale_factor, residuals((length(locs_InSAR)+2):3:end)'*scale_factor, 'k', 'Autoscale', 'off');   % i.e. plot every third row of the locations against E and N...
        end



        axis equal
        colorbar

        % plot fault trace
        for j = 1 : n_fault_strands
                x = [fault_coords(j,1), fault_coords(j,3)];
                y = [fault_coords(j,2), fault_coords(j,4)];
                hold on;
                plot(x, y, 'm', 'Linewidth', 2)
        end

        %xlabel('UTM x, I think')
        %ylabel('UTM y, I think')
        if counter ==3
            title('Residuals');
        end
        

        [redbluecmap] = redblue;
        colormap(flipud(redbluecmap))


        % Make sure all colorbars are same scale, only display on RHS
        if strcmp(data.InSAR_datafile, 'none') ~= 1
            MaxC = max( abs([max(d_InSAR), max(d_hat(1:length(locs_InSAR))), max(residuals(1:length(locs_InSAR)))]));  % max colour is max InSAR value - ignore GPS arrows
            %MinC = min( [min(d_InSAR), min(d_hat(1:length(locs_InSAR))), min(residuals(1:length(locs_InSAR)))]);  % max colour is max InSAR value - ignore GPS arrows
            %subplot(max([n_InSAR_scenes 1]),3,counter-2); caxis([ MinC MaxC]); colorbar('off');
            %subplot(max([n_InSAR_scenes 1]),3,counter-1); caxis([ MinC MaxC]); colorbar('off');
            %subplot(max([n_InSAR_scenes 1]),3,counter); caxis([ MinC MaxC]);
            subplot(max([n_InSAR_scenes 1]),3,counter-2); caxis([ -MaxC MaxC]); colorbar('off');
            subplot(max([n_InSAR_scenes 1]),3,counter-1); caxis([ -MaxC MaxC]); colorbar('off');
            subplot(max([n_InSAR_scenes 1]),3,counter); caxis([ -MaxC MaxC]);
            %subplot(3,1,1);  caxis([ MinC MaxC]);
            %subplot(3,1,2);  caxis([ MinC MaxC]); colorbar('off');
            %subplot(3,1,3);  caxis([ MinC MaxC]); colorbar('off');
        end
        if strcmp(data.GPS_datafile, 'none') ~= 1
            ylabel(colorbar, 'LOS displacement (m)')    % note that GPS has been projected into LOS
        end
        colormap(flipud(redbluecmap))
        if strcmp(data.atolls_datafile, 'none') ~= 1
            ylabel(colorbar, 'vertical displacement (m)')    % note that GPS has been projected into LOS
        end

        %linkaxes([ax1, ax2, ax3], 'xy')

        counter = counter+1;




      end


end




%% Plot histograms for all slip patches

if strcmp( display.plothists, 'plothistsall') == 1
    % Subplot numbers the plots along the rows. But I have numbered my
    % patches down the columns. So I have to create a matrix re-numbering
    % my slip patches in this way. So that:
    %
    %  Patch numbers:                  Subplot numbers:
    %   1 3 5               becomes     1 2 3
    %   2 4 6                           4 5 6
    %
    % So that patch 1 is subplot 1, but patch 2 is subplot 4, etc

    A = 1: total_n_slip_patches;
    A = reshape(A, total_n_along_strike_patches, total_n_down_dip_patches);
    A = A';
    A = fliplr(A);
    plot_numbers = reshape(A, total_n_slip_patches, 1);

    nbins = 20;
    
    % Plot the patches
    figure('position', [100, 350, 1600, 1200])
    for i = 1: total_n_slip_patches
       subplot( total_n_down_dip_patches, total_n_along_strike_patches, plot_numbers(i))
       colormap default
       %hist( slip_keep_just_real_patches(i,:), nbins)
       %[counts, centers] = hist(slip_keep_just_real_patches(i,:), nbins);   
       hist(slip_keep(i,onlyonkeptpatches(i,:)), nbins);
       hold on
       x_mean = [patch_mean(i) patch_mean(i)];
       x_mode = [patch_mode(i) patch_mode(i)];
       x_median = [patch_median(i) patch_median(i)];
       %x_conf_lower = [patch_conf_intervals(i,1) patch_conf_intervals(i,1)];
       %x_conf_higher = [patch_conf_intervals(i,2) patch_conf_intervals(i,2)];
       plot(x_mean,y,'r');
       plot(x_mode,y,'g');
       plot(x_median,y,'k');
       %plot(x_conf_lower, y, 'y');
       %plot(x_conf_higher, y, 'y');
           if strcmp(testing.testing_mode, 'yes')
              load(testing.making_model, 'synthetic_slip');
              if length(synthetic_slip) == total_n_slip_patches;
                x_true = [synthetic_slip(i) synthetic_slip(i)]; % plot x_true as a line
                plot(x_true, y, 'm');
              end
           end
       set(gca,'YTick',[]);
       set(gca,'Xticklabel',[])
       %title(['Slip patch ', num2str(i)])
       ylim([0 invert.iterations/3])
       xlim([0 max(max(slip_keep))])
    end
  
% Put x and y labels on middle-outside plots
subplot(total_n_down_dip_patches, total_n_along_strike_patches, plot_numbers(ceil(total_n_down_dip_patches/2)))
ylabel('Frequency')
subplot(total_n_down_dip_patches, total_n_along_strike_patches, plot_numbers((ceil(n_along_strike_patches/2))*total_n_down_dip_patches))
xlabel('Slip (m)')
          

    % Put a legend, but only on the histogram on the far right
    subplot(total_n_down_dip_patches, total_n_along_strike_patches, total_n_along_strike_patches)
    if strcmp(testing.testing_mode, 'yes')
        legend('Sampling', 'Mean', 'Mode', 'Median','95%','95%','True');
    else
        legend('Sampling', 'Mean', 'Mode', 'Median','95%','95%');
    end
    legend('location', 'northeastoutside')
end


%% Plot pdfs of correlation lengths, if necessary

if strcmp(invert.smoothing, 'VK') ==1 && strcmp(invert.solve_for_correlation_length, 'yes') == 1
   
    figure('position', [100, 300, 600, 800])
    subplot(2,1,1)
    hist(a_as_keep, 100)
    title('model parameter a along-strike')
    ylabel('frequency')
    
    subplot(2,1,2)
    hist(a_dd_keep, 100)
    title('model parameter a down-dip')
    xlabel('value')
    
end


%% Plot pdfs of dips, if necessary

if strcmp(invert.solve_for_dip, 'yes') == 1;
    
    figure('position', [100, 300, 600, 800])
    title('dip')
    
    if strcmp(testing.testing_mode, 'yes') == 1;
        true_model = load(testing.making_model, 'synthetic_disloc_model');
        dip_true = true_model.synthetic_disloc_model(4,:);
    end
    
    clear counts
    clear centers
    
    for i = 1: n_fault_strands
        


        
%          nbins2 = 100;
%          for r = 1: total_n_slip_patches
%              [counts2(r,:), centers2(r,:)] = hist(dip_keep(i,r), nbins2);   % count the number in each bin and the center of each bin
%          end
%          [M,I] = max(counts2');
    
       subplot(n_fault_strands, 1, i);
       hist(dip_keep(i,:), 100);
       title(['dip on fault plane ', num2str(i)])
       ylabel('frequency')
       
       
            % Find mode

                [counts, centers] = hist(dip_keep(i,:), 100);   % count the number in each bin and the center of each bin
                [M,I] = max(counts'); % If A is a matrix, then max(A) is a row vector containing the maximum value of each column
                
                dip_mode = centers(1,I);

                hold on;
                x = [dip_mode dip_mode];
                y = [0 max(counts)];
                plot( x,y, 'r', 'Linewidth', 5)
    
       
       if strcmp(testing.testing_mode, 'yes') == 1;
            x = [dip_true(first_patch_in_strand(i)) dip_true(first_patch_in_strand(i))];
            y = [0 500]; %y = [0 max(M)];
            hold on;
            plot(x,y,'r', 'Linewidth', 2);
       end
        
        if i == n_fault_strands;
            xlabel('dip')
        end
        
    end

  clear counts
  clear centers   
    
end

    




%% Plot in 3D

if strcmp(display.plot3d, 'yes') == 1

    h=figure('position', [200, 300, 1600, 1000]);
    set(gca,'fontsize',30)
    
    % PROBLEM: gareth funning's doplot3d needs utmx utmy to calculate
    % sizes. But I have it in lat long and it's in different utm zones...
    
    % Plot fault in 3D with MAP solution
    faults = disloc_model;
    
    % translate back to lat and long
        if strcmp(use_local_coordinate_system, 'yes') == 1
            xy(1,:) = disloc_model(1,:)/1000; 
            xy(2,:) = disloc_model(2,:)/1000; 
            [longlat] = local2llh(xy,data.origin);  % this is fine, it's just being translated back 
            
            % translate to utmx utmy
            [utmx, utmy, UTMzone] = ll2utm(longlat(2,:), longlat(1,:));
            faults(1,:) = utmx;
            faults(2,:) = utmy;
        end
    
    faults(5,:) = rake_mostlikely;
    faults(6,:) = patch_mean;
    %faults(6,:) = patch_std;
    doplot3d(faults', 'hot');
    %doplot3d_ruthhack_utm2ll(faults', 'jet', data.UTMzone*ones(total_n_slip_patches,1))
    hold on;
    colorbar
    %caxis([0, cmaxvalue]);
    caxis([0, max(faults(6,:))]);
    caxis([0, 3]);
    title('3D mean')
    xlabel('UTM x')
    ylabel('UTM y')
    zlabel('Depth (km)')
    ylabel(colorbar, 'Slip (m)', 'FontSize', 30)
    hot2 = hot;
    hot2(1:40,:) = [];      % chopping off start of colorbar to get rid of black - i want the maximum colour to be dark red  
    hot2 = [hot2; ones(30,3)];  % adding more white to the bottom
    colormap(hot2)
    colormap(flipud(colormap))
    
    %h = get(fig,'CurrentAxes');
    %direction = [ 1 0 0 ];
    %rotate(h, direction, 25);
    
    % Plot epicenter
    %[utmx, utmy] = ll2utm(data.EQ_epicenter(1), data.EQ_epicenter(2));
    %scatter3( utmx, utmy, data_EQ_epicenter(3), 'yp', 'filled');
    
    % Add rake
    quiv_mags = [faults(6,:)'.*cosd(rake_mode), faults(6,:)'.*sind(rake_mode)];
    [xquivmag,yquivmag,zquivmag] = reproject_quiv(quiv_mags,disloc_model(3,:)',disloc_model(4,:)');   % nicked from Tom.   reproject_quiv(ss_mag, ds_mag, strike, dip)
    %[xquivmag,yquivmag,zquivmag] = pol2cart(disloc_model(3,:)',rake_mode,patch_mode);
    scale_factor = 0.7;
    
    % work out where to plot rake points
    dip = disloc_model(4,1);
    extrawidth = abs(cosd(dip) * (disloc_model(8,:)+disloc_model(9,:))'/(-2*sind(disloc_model(4,1))))';
    x_extra = abs(extrawidth * cosd(strike));
    y_extra = abs(extrawidth * sind(strike));
    hold on;
    if strike > 90 && strike <= 180 || strike > 180 && strike <= 270
         x_extra = -x_extra;
    end
    if strike >= 0 && strike <= 90 || strike > 90 && strike <= 180
        y_extra = -y_extra;
    end
    %quiver3(   (disloc_model(1,:)'+ x_extra )/1000,     (disloc_model(2,:)'- y_extra)/1000 ,    (disloc_model(8,:)+disloc_model(9,:))'/(-2*sind(disloc_model(4,1)))/1000    ,xquivmag*scale_factor, yquivmag*scale_factor, zquivmag*scale_factor, 'k', 'Linewidth', 1.5, 'Autoscale', 'off');
    quiver3(   (disloc_model(1,:) + x_extra )/1000,     (disloc_model(2,:) + y_extra)/1000 ,    -((disloc_model(9,:)-disloc_model(8,:))/2 + disloc_model(8,:))/1000,       xquivmag'*scale_factor,      yquivmag'*scale_factor,       zquivmag'*scale_factor, 'k', 'Linewidth', 1.5, 'Autoscale', 'off');
   
end

         

%% Plot 2D pdfs
% With many thanks to David Bekaert for this script

if strcmp(display.plotmarginalPDFs, 'yes') == 1
 
        slip_for_marginal_PDF_plotting = slip_keep;
    
                %true_slip_for_marginal_PDF_plotting = synthetic_slip;
    
     if total_n_slip_patches > 15      % if you have a lot of slip patches, plot the patches with the highest slip only   
           %disp('JUST A NOTE - since you have so many slip patches, your marginal PDFs plot will look silly...');       
%          random_selection = round((total_n_slip_patches).*rand(1000,1));  % pick a few slip patches randomly to plot against each other
%             while any(random_selection == 0)
%                 random_selection = round((total_n_slip_patches).*rand(15,1));
%             end
%          for k = random_selection;
%             slip_for_marginal_PDF_plotting = slip_for_marginal_PDF_plotting(k,:);
%          end 
            [~,sortIndex] = sort(patch_MAP,'descend');
            maxIndex = sortIndex(1:5);
            patches_with_max_slip = slip_for_marginal_PDF_plotting(maxIndex, :);
            slip_for_marginal_PDF_plotting = patches_with_max_slip;
           if strcmp(testing.testing_mode, 'yes') == 1
                true_slip_for_marginal_PDF_plotting = synthetic_slip(maxIndex);
           end
     else
         sortIndex = 1:total_n_slip_patches;
        if strcmp(testing.testing_mode, 'yes') == 1
                true_slip_for_marginal_PDF_plotting = synthetic_slip;
        end
     end
           
     
    if size(slip_for_marginal_PDF_plotting, 2) > size(slip_for_marginal_PDF_plotting, 1)   % transpose if necessary
        slip_for_marginal_PDF_plotting = slip_for_marginal_PDF_plotting';
    end
    
    
    figure('position', [100, 300, 1600, 1200])
    if strcmp(testing.testing_mode, 'yes') == 1
        [H,AX,BigAx,P,PAx,cc_colorbar] = plotmatrix_lower2(slip_for_marginal_PDF_plotting,'plot_color', true_slip_for_marginal_PDF_plotting);
    else
         [H,AX,BigAx,P,PAx,cc_colorbar] = plotmatrix_lower(slip_for_marginal_PDF_plotting,'plot_color');        % with many thanks to David Bekaert for this script %   [H,AX,BigAx,P,PAx,cc_colorbar] = plotmatrix_lower_david(data(ix_range_3,:),'plot_color',data_opt);   
    end
    

    for i = 1: size(slip_for_marginal_PDF_plotting,2)
       strings{i} = (['P', num2str(sortIndex(i))]) ;
       strings_units{i} = (['P', num2str(sortIndex(i)), ' slip (m) ']);                      % (['Patch ', num2str(i), '  ';'slip [m] '])
    end


    for k = 1 : (size(slip_for_marginal_PDF_plotting,2)-1)                      % loop to (n-1) because we want to put the labels on the LHS from slip patch 2:end, not 1:end. Then add on the bottom label for last slip patch.
       set(get(PAx(k),'xlabel'),'string',strings_units{k},'fontsize',fontsize_plot-6);         %set(get(PAx(k),'xlabel'),'string',strings_units{k},'fontsize',fontsize_plot-1)
       set(get(AX(k,1),'ylabel'),'string',strings{k+1},'fontsize',fontsize_plot-6);     % set(get(AX(k,1),'ylabel'),'string',strings{k+1},'fontsize',fontsize_plot-1)
    end
    set(get(PAx( size(slip_for_marginal_PDF_plotting,2)),'xlabel'),'string',strings_units{size(slip_for_marginal_PDF_plotting,2)},'fontsize',fontsize_plot-6);      % add on bottom label for last slip patch, name is the last value in the matrix 'strings'
    %set(get(AX(1), 'title'), 'string', 'Patch 1','fontsize',fontsize_plot-1);
    %set(get(BigAx, 'title'), 'string', 'Marginal PDFs','fontsize',fontsize_plot-1);

%     
%     if strcmp(testing.testing_mode, 'yes') == 1;
%         
%         plot_number_to_change = [1 6 7 11 12 13 16 17 18 19];
%         
%         for i = 1:5
%             for j = 1:5;
%                 subplot(5, 5, plot_number_to_change(i+j));
%                 scatter(true_slip_for_marginal_PDF_plotting(i), true_slip_for_marginal_PDF_plotting(j));
%             end
%             
%         end
%         
%     end
    
    

end


%%  Plot parameters through time


% if strcmp(invert.inversion_type, 'bayesian') == 1
% 
%     % Slip
%     figure('position', [200, 300, 1600, 1200])
%     for i = 1: total_n_slip_patches
%        subplot(total_n_down_dip_patches,total_n_along_strike_patches,plot_numbers(i))
%        plot(1:length(slip_keep), slip_keep(i,:));
%        ylim([priors.min_slip(1),max(max(slip_keep))])
%     end
%     xlabel('Iterations')
%     subplot(total_n_down_dip_patches,total_n_along_strike_patches,1)
%     ylabel('Slip')
%     title('Slip trace plot')
% 
%     if strcmp(invert.smoothing, 'none') ~= 1
%         % Alpha
%         figure('position', [200, 300, 800, 1200])
%         for i = 1: n_fault_strands_for_smoothing
%            subplot(n_fault_strands_for_smoothing,1,i)
%            plot(1:length(alpha2_keep), alpha2_keep(i,:));
%            ylim([priors.min_alpha2(i),max(alpha2_keep(i,:))])
%         end
%         xlabel('Iterations')
%         subplot(n_fault_strands_for_smoothing,1,1);
%         ylabel('Alpha')
%         title('Alpha trace plot')
%     end
% 
%     % Rake
%     figure('position', [200, 300, 1600, 1200])
%     for i = 1: total_n_slip_patches
%        subplot(total_n_down_dip_patches,total_n_along_strike_patches,plot_numbers(i))
%        plot(1:length(rake_keep), rake_keep(i,:));
%        ylim([priors.min_rake(1),priors.max_rake(1)])
%     end
%     xlabel('Iterations')
%     subplot(total_n_down_dip_patches,total_n_along_strike_patches,1)
%     ylabel('Rake')
%     title('Rake trace plot')
% 
%     % Circharm parameters
% 
% end

%% Offset

if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
   
    figure('position', [600, 300, 800, 1200])
    title('Offset')
    for p = 1: n_InSAR_scenes
        subplot(n_InSAR_scenes,1,p)
        hist(offset_keep(p,:), 50)
    end
end
 

if strcmp(invert.solve_for_InSAR_ramp, 'yes') == 1
    figure('position', [600, 300, 800, 1200])
    for p = 1: n_InSAR_scenes
        subplot(n_InSAR_scenes,3,3*p -2)
        hist(ramp_keep(3*p -2,:), 50)
        ylabel(['InSAR scene ', num2str(p)]);
        if p == 1
            title('InSAR ramp x')
        end

        subplot(n_InSAR_scenes,3,3*p -1)
        hist(ramp_keep(3*p -1,:), 50)
        if p == 1
            title('InSAR ramp y')
        end

        subplot(n_InSAR_scenes,3,3*p)
        hist(ramp_keep(3*p,:), 50)
        if p == 1
            title('InSAR ramp offset')
        end 
    end
end

%%  Plot the walkers


if strcmp(invert.ensemble_sampling, 'yes') == 1
    
   figure('position', [200, 300, 1600, 1200]) 
    
%    % Slip
%    subplot(2,2,1)
%    randomslippatch = ceil(total_n_slip_patches*rand);
%    hold all
%    plot(alpha2walkers_keep, 1:length(alpha2walkers_keep));
%    xlabel(['slip patch ', num2str(randomslippatch)])
        figure
        randomslippatches = ceil(total_n_slip_patches*rand(6,1));
        hold all
        slipwalkers = m_keep(m_identifyer==1,:,:);
        n = 1;
        slipplotwalkers = zeros(size(m_keep,2),nwalkers);
        for p = 1:6
            subplot(2,3,n)
            slipplotwalkers(:,:) = slipwalkers(randomslippatches(p), :,:);
            plot(slipplotwalkers, 1:length(slipplotwalkers));
            n = n+1;
        end
        xlabel('slip')
        ylabel('i')
    
%    % Rake
        figure
        randomrakepatches = ceil(total_n_slip_patches*rand(6,1));
        hold all
        rakewalkers = m_keep(m_identifyer==3,:,:);
        n = 1;
        rakeplotwalkers = zeros(size(m_keep,2),nwalkers);
        for p = 1:6
            subplot(2,3,n)
            rakeplotwalkers(:,:) = rakewalkers(randomrakepatches(p), :,:);
            plot(rakeplotwalkers, 1:length(rakeplotwalkers));
            n = n+1;
        end
        xlabel('rake')
        ylabel('i')
   
   % Alpha
   figure
   hold all
   alpha2walkers = zeros(size(m_keep,2),nwalkers);
   alpha2walkers(:,:) = m_keep(m_identifyer==2,:,:);
   plot(alpha2walkers, 1:length(alpha2walkers));
   xlabel('alpha^2')
   ylabel('i')
   
   % Offset
   if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
       figure
       title('InSAR offset')
       n=1;
       hold all
       offsetwalkers =  zeros(size(m_keep,2),nwalkers);
       offsetwalkers(:,:) = m_keep(m_identifyer==2,:,:);
       for p= 1:n_InSAR_scenes
           subplot(n_InSAR_scenes, 1, n)
           plot(alpha2walkers, 1:length(alpha2walkers));
           n=n+1;
       end
       xlabel('offset')
       ylabel('i')
   end

end

