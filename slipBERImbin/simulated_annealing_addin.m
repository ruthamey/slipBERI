        % Initial values
        m_simulatedannealing = m_initial(m_identifyer ~= 7 & m_identifyer ~= 8 & m_identifyer ~= 9 & m_identifyer ~= 10 & m_identifyer ~= 11);  % not trying to solve for circ harmonic parameters
        step_sizes_simulatedannealing = step_sizes(m_identifyer ~= 7 & m_identifyer ~= 8 & m_identifyer ~= 9 & m_identifyer ~= 10 & m_identifyer ~= 11);
        m_min_simulatedannealing = m_min(m_identifyer ~= 7 & m_identifyer ~= 8 & m_identifyer ~= 9 & m_identifyer ~= 10 & m_identifyer ~= 11);
        m_max_simulatedannealing = m_max(m_identifyer ~= 7 & m_identifyer ~= 8 & m_identifyer ~= 9 & m_identifyer ~= 10 & m_identifyer ~= 11);
        T_trials = [1000 100 10 1 0.1 0.01 0.001 0.0001];
        isim = 10000;
        rejections = zeros(1, length(T_trials));
        m_keep = zeros(length(m_simulatedannealing),isim);
        sigma_keep = zeros(isim, length(T_trials));
        
        for n = 1: length(T_trials)
            
            T = T_trials(n);
            
            G_curr = G_ss*diag(cosd(m_simulatedannealing(m_identifyer==3,n))) +  G_ds*diag(sind(m_simulatedannealing(m_identifyer==3,n)));  
           
            % Sort out alpha^2
            if strcmp(priors.alpha2_prior, 'logarithmic') == 1
                alpha2_curr = 10.^(m_curr(m_identifyer==2));
            else
                alpha2_curr = m_curr(m_identifyer==2);
            end
            
            % Sort out offset
            if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
                offset_modelparameter_temp = m_trial(m_identifyer==5);
                for q = 1:n_InSAR_scenes
                    offset_temp(InSAR_identifyer==q,1) = offset_modelparameter_temp(q);
                end
                offset_temp((end+1):length(d))=0;    % No offset on GPS
            else
                offset_temp = zeros( length(d), 1);
            end
            
            [prior_curr, ~]=  calc_logprior_VK( m_curr(m_identifyer==1), inv_sigma_s, alpha2_curr, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing);
            likelyhood_curr= calc_loglikely(m_curr(m_identifyer==1), d, G_curr, inv_sigma_d, offset_temp, m_curr(m_identifyer==6) );
            %L_unscaled = (2*pi*alpha_curr)^(length(m_simulatedannealing)/2) * exp(prior+likelyhood) ;
            %L = (2*pi*alpha_curr)^(length(m_simulatedannealing)/2) * exp((prior+likelyhood)/ T) ;
            logL_curr_unscaled = log((2*pi*alpha2_curr)^(-length(m_simulatedannealing(:,n))/2)) + ((prior_curr+likelyhood_curr)) ;
            
            for i = 1:isim
                
                % Add step size
                m_trial = m_simulatedannealing(:,n) + ((step_sizes_simulatedannealing - (-step_sizes_simulatedannealing)) .*rand(length(m_simulatedannealing(:,n)),1) - step_sizes_simulatedannealing);
                
                % Check in priors
                toosmall = m_trial < m_min_simulatedannealing;
                m_trial(toosmall) = 2*m_min_simulatedannealing(toosmall) - m_trial(toosmall);
                toobig = m_trial > m_max_simulatedannealing;
                m_trial(toobig) =  2*m_max_simulatedannealing(toobig) - m_trial(toobig);
                
                G_trial = G_ss*diag(cosd(m_trial(m_identifyer==3))) +  G_ds*diag(sind(m_trial(m_identifyer==3)));  
              
                % Sort out alpha^2
                if strcmp(priors.alpha2_prior, 'logarithmic') == 1
                    alpha2_trial = 10.^(m_trial(m_identifyer==2));
                else
                    alpha2_trial = m_trial(m_identifyer==2);
                end
                
                % Sort out offset
                if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
                    offset_modelparameter_temp = m_trial(m_identifyer==5);
                    for q = 1:n_InSAR_scenes
                        offset_temp(InSAR_identifyer==q,1) = offset_modelparameter_temp(q);
                    end
                    offset_temp((end+1):length(d))=0;    % No offset on GPS
                else
                    offset_temp = zeros( length(d), 1);
                end
                
                [prior_trial, ~]= calc_logprior_VK( m_trial(m_identifyer==1), inv_sigma_s, alpha2_trial, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing);
                likelyhood_trial= calc_loglikely(m_trial(m_identifyer==1), d, G_trial, inv_sigma_d, offset_temp, 1 );
                %L_trial = (2*pi*alpha_trial)^(length(m_simulatedannealing)/2) * exp( (prior+likelyhood)/ T);
                logL_trial_unscaled = log((2*pi*alpha2_trial)^(-length(m_simulatedannealing(:,n))/2)) + (prior_trial+likelyhood_trial);
                
                
                prior_ratio = exp( (prior_trial - prior_curr));
                alpha2_ratio = prod( (alpha2_trial ./ alpha2_curr) .^ (-n_slip_patches_ON_on_each_fault_strand_for_smoothing/2));      % do the product because we need to multiply all the ratios for the n fault strands
                likelihood_ratio = exp( (likelyhood_trial - likelyhood_curr)/T);
                ratio = prior_ratio *alpha2_ratio * likelihood_ratio ;
                
                % If L_trial > L, or L_trial/ L > random number, keep new values.
                if ratio > rand
                    %L = L_trial;                        % Update L
                    m_simulatedannealing(:,n) = m_trial;     % Update m_simulatedannealing
                    prior_curr = prior_trial;
                    likelyhood_curr = likelyhood_trial;
                    alpha2_curr = alpha2_trial;
                    logL_curr_unscaled = logL_trial_unscaled;
                else                                    % m_simulatedannealing, L stay the same
                    rejections(n) = rejections(n) + 1;  % This is a rejection. Shot down.
                end
                
                m_keep(:,i) = m_simulatedannealing(:,n);     % Store trial value of m_simulatedannealing (either updated version or old version, depending on whether it passed the Metropolis rule.)
                sigma_keep(i,n) = logL_curr_unscaled;  % store unscaled value of sigma
                
            end
            
            % Display the number of rejections
            disp(['Number of rejections for T = ',num2str(T),' is ', num2str(rejections(n)),])
            
            % Find m_simulatedannealing values with max a posteriori to use as starting values for next iteration
            [~, pos] = max(sigma_keep(:,n));
            m_simulatedannealing(:,n+1) = m_keep(:,pos);     % This is assigned to m in the start of the next loop
            
        end
        
        % Make initial values
        m_initial(m_identifyer ~= 7 & m_identifyer ~= 8 & m_identifyer ~= 9 & m_identifyer ~= 10 & m_identifyer ~= 11) = m_simulatedannealing(:,end);% the 'best' values are the ones from the last simulated annealing iteration
        m_curr = m_initial;
        m_trial = m_curr;
        if strcmp(priors.alpha2_prior, 'logarithmic') == 1
            alpha2_curr = 10.^(m_curr(m_identifyer==2));
            alpha2_initial = alpha2_curr;
        else
            alpha2_curr = m_curr(m_identifyer==2);
            alpha2_initial = alpha2_curr;
        end
        if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
            offset_modelparameter_temp = m_trial(m_identifyer==5);
            for q = 1:n_InSAR_scenes
                offset_curr(InSAR_identifyer==q,1) = offset_modelparameter_temp(q);
            end
            offset_curr((end+1):length(d))=0;    % No offset on GPS
        else
            offset_curr = zeros( length(d), 1);
        end
        
        
% %********** Solve for which circular harmonic variables *******************
%         if strcmp(invert.circular_harmonics, 'yes') == 1
%             
%             disp('  ')
%             disp('Solving for initial circharm parameters using simulated annealing result.')
%             disp('  ')
%             
%             % Turn off patches with less than 1/3 of maximum slip
%             onoffidentifyer_aim = onoffidentifyer;
%             onoffidentifyer_aim(m_initial(m_identifyer==1)<(max(m_initial(m_identifyer==1))/3)) = 0;   % Turn off patches with less than 1/3 max slip
%               
%             % Solve for circharm parameters that can turn off these patches
%             icircharm = 10000;
%             ch_coeffs_curr = m_initial(m_identifyer ==7);
%             ch_phi_curr = m_initial(m_identifyer ==8);
%             ch_center_curr = m_initial(m_identifyer ==9);
%             ch_coeffs_keep = zeros(n_harmonics,icircharm);
%             ch_phi_keep = zeros(n_harmonics-1,icircharm);
%             ch_center_keep = zeros(2,icircharm);
%             logprob_keep = zeros(1,icircharm);
%             %ch_sigma_d = diag(((max(m_initial(m_identifyer==1)))./m_initial(m_identifyer==1)));  % Give more weight to squares that have higher slip (we REALLY want these patches to be on)
%             ch_sigma_d = diag( ones(length(m_initial(m_identifyer==1)),1));
%             inv_ch_sigma_d = inv(ch_sigma_d);
%             logprob_curr = -(onoffidentifyer_aim - onoffidentifyer)'*inv_ch_sigma_d*(onoffidentifyer_aim - onoffidentifyer);
%             
%             for qq = 1:icircharm
%                 
%                 % Add random number to circharm parameters
%                 ch_coeffs_trial = ch_coeffs_curr + (-250 + (500).*rand(n_harmonics,1));      % plus or minus 250
%                 ch_phi_trial = ch_phi_curr + (-0.5 + (1).*rand(n_harmonics-1,1));            % plus or minus 1 radian
%                 ch_center_trial = ch_center_curr + (-250 + (500).*rand(2,1));               % plus or minus 250
%                 
%                 % Make sure not too small / too big   - THIS IS CLUNKY
%                 toosmall = ch_coeffs_trial < min_circharm_coeffs;
%                 ch_coeffs_trial(toosmall) = 2*min_circharm_coeffs(toosmall) - ch_coeffs_trial(toosmall);
%                 toobig = ch_coeffs_trial > max_circharm_coeffs;
%                 ch_coeffs_trial(toobig) =  2*max_circharm_coeffs(toobig) - ch_coeffs_trial(toobig);
%                 %
%                 toosmall = ch_phi_trial < min_circharm_phi;
%                 ch_phi_trial(toosmall) = 2*min_circharm_phi(toosmall) - ch_phi_trial(toosmall);
%                 toobig = ch_phi_trial > max_circharm_phi;
%                 ch_phi_trial(toobig) =  2*max_circharm_phi(toobig) - ch_phi_trial(toobig);
%                 %
%                 toosmall = ch_center_trial < min_circharm_center;
%                 ch_center_trial(toosmall) = 2*min_circharm_center(toosmall) - ch_center_trial(toosmall);
%                 toobig = ch_center_trial > max_circharm_center;
%                 ch_center_trial(toobig) =  2*max_circharm_center(toobig) - ch_center_trial(toobig);
%                 
%                 % Work out which patches are 'on' and 'off'
%                 [circx,circz]=circharm(ch_coeffs_trial,ch_phi_trial,0);     % Set to '1' if want to plot
%                 circharm_center_trial = ch_center_trial;
%                 onoffidentifyer_trial = inpolygon(patchx, patchz, (circx+circharm_center_trial(1)), (circz+circharm_center_trial(2)))';
%                 
%                 % Calculate new probability
%                 logprob_trial = -(onoffidentifyer_aim - onoffidentifyer_trial)'*inv_ch_sigma_d*(onoffidentifyer_aim - onoffidentifyer_trial);                % log (   exp (d-Gm)'(d-Gm) )
%                 
%                 if exp( logprob_trial - logprob_curr) > rand
%                     ch_coeffs_curr = ch_coeffs_trial;
%                     ch_phi_curr = ch_phi_trial;
%                     ch_center_curr = ch_center_trial;
%                     logprob_curr = logprob_trial;
%                     
%                 end
%                 
%                 ch_coeffs_keep(:,qq) = ch_coeffs_curr;
%                 ch_phi_keep(:,qq) = ch_phi_curr;
%                 ch_center_keep(:,qq) = ch_center_curr;
%                 logprob_keep(:,qq) = logprob_curr;
%             end
%             
%             % Find most likely solution
%             [~, I] = max(logprob_keep);
%             circharm_coeffs_initial = ch_coeffs_keep(:, I);
%             circharm_phi_initial = ch_phi_keep(:,I);
%             circharm_center_initial = ch_center_keep(:,I);
%             
% %             % start with circharm bigger than whole fault
% %             circharm_center_initial = [fault_length_for_smoothing/2; fault_width_for_smoothing/2];       % x,z offset
% %             circharm_coeffs_initial = [ sqrt((fault_length_for_smoothing/2)^2 + (fault_width_for_smoothing/2)^2); zeros(n_harmonics-1,1)];      % so the circle starts bigger than the fault
% %             circharm_phi_initial = zeros(n_harmonics-1,1);
%             
%             % start with right-ish area
%             %     circharm_center_initial = [fault_length_for_smoothing/2; fault_width_for_smoothing/2];       % x,z offset
%             %     circharm_coeffs_initial = [3000; 0; 3600; 0];
%             %     circharm_phi_initial = [0; (pi/3+pi/2); 0];
% 
%             % Check
%             [circx,circz]=circharm(circharm_coeffs_initial,circharm_phi_initial,0);     % Set to '1' if want to plot
%             onoffidentifyer = inpolygon(patchx, patchz, (circx+circharm_center_initial(1)), (circz+circharm_center_initial(2)))';
%             
%             figure;
%             scatter(patchx(onoffidentifyer_aim==1), patchz(onoffidentifyer_aim==1), 1800, 'gs', 'filled')   % plot x and z of on patches
%             hold on;
%             scatter(patchx(onoffidentifyer_aim==0), patchz(onoffidentifyer_aim==0),1800, 'rs', 'filled')   % plot x and z of off patches
%             plot(circx+circharm_center_initial(1),circz+circharm_center_initial(2));
%             set(gca,'Ydir','reverse')
%             axis equal tight
%             axis([min([circz+circharm_center_initial(1) 0]) max([circz+circharm_center_initial(1)  fault_length_for_smoothing]) min([circx+circharm_center_initial(2) 0]) max([circx+circharm_center_initial(2)  fault_width_for_smoothing])])
%             ylabel('Down dip')
%             xlabel('Along strike')
%             title('Slipping area on fault')
%             
%             % Use these values as starting values - recalcualte everything
%             m_initial(m_identifyer_master == 7) = circharm_coeffs_initial;
%             m_initial(m_identifyer_master == 8) = circharm_phi_initial;
%             m_initial( m_identifyer_master == 9) = circharm_center_initial;
%             mcircharmfreeze = m_initial(circharmparameters==1);
%             
%             for n = 1:n_fault_strands_for_smoothing
%                 n_slip_patches_ON_on_each_fault_strand_for_smoothing(n,1) = sum(onoffidentifyer(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),1));
%             end
%             previous_n_slip_patches_ON_on_each_fault_strand_for_smoothing = n_slip_patches_ON_on_each_fault_strand_for_smoothing;
%             first_patch_in_strand_for_smoothing(1) = 0;
%             first_patch_in_strand_for_smoothing(2:(n_fault_strands_for_smoothing+1),1) = cumsum(n_slip_patches_ON_on_each_fault_strand_for_smoothing);
%             first_patch_in_strand_for_smoothing = first_patch_in_strand_for_smoothing+1;
%             first_patch_in_strand_for_smoothing(n_fault_strands_for_smoothing+1) = [];
%             last_patch_in_strand_for_smoothing = cumsum(n_slip_patches_ON_on_each_fault_strand_for_smoothing);
%             
%             %  Update fault size - NOTE THIS ONLY WORKS FOR ONE FAULT STRAND FOR NOW. If doing for more than one fault strands need to sum along the fault, not just the change in utmx and utmy between the patches.
%             tippytop = min(disloc_model(8,onoffidentifyer==1));     % this finds the minimum top depth of any patches that are on from the matrix disloc model
%             verybottom = max(disloc_model(8,onoffidentifyer==1));   % this finds the maximum top depth of any patches that are on from the matrix disloc model
%             slippingpatch_width_for_smoothing =  verybottom-tippytop;
%             
%             farleft = min(disloc_model(1,onoffidentifyer==1));
%             farright = max(disloc_model(1,onoffidentifyer==1));
%             slippingpatch_length_for_smoothing = abs(farleft-farright);           % abs in case in a synthetic test you've defined +x and -x
%             
%             if length(slippingpatch_width_for_smoothing) == 0
%                 slippingpatch_width_for_smoothing = 0;
%             end
%             
%             if length(slippingpatch_length_for_smoothing) == 0
%                 slippingpatch_length_for_smoothing = 0;
%             end
%             
%             % Recalculate correlation lengths - drawing from distribution if solving for correlation lengths.            
%             if strcmp(invert.solve_for_correlation_length, 'yes') == 1
%                 if strcmp(predominant_faulting_style, 'ss') == 1
%                     % Strike slip parameters
%                     a_as_ss_mean = 1860 + 0.34*slippingpatch_length_for_smoothing;
%                     a_as_ss_std = sqrt(1120^2+0.03^2*slippingpatch_length_for_smoothing^2);          % Errors from (Mai and Beroza, 2002), Table 2, Pg 14
%                     a_dd_ss_mean = -390 + 0.44 *slippingpatch_width_for_smoothing;                   % Propagation of errors from talbe on Wikipedia (shh)  https://en.wikipedia.org/wiki/Propagation_of_uncertainty
%                     a_dd_ss_std = sqrt(470^2+0.04^2*slippingpatch_width_for_smoothing^2);
%                 elseif strcmp(predominant_faulting_style, 'ds') == 1
%                     % Dip slip parameters
%                     a_as_ds_mean =  1100 + 0.31*slippingpatch_length_for_smoothing;
%                     a_as_ds_std = sqrt(400^2+0.01^2*slippingpatch_length_for_smoothing^2);
%                     a_dd_ds_mean =  580 + 0.35*slippingpatch_width_for_smoothing;
%                     a_dd_ds_std = sqrt(240^2+0.01^2*slippingpatch_width_for_smoothing^2);
%                 end
%             end
%             
%             for n = 1 : n_fault_strands_for_smoothing
%                 if strcmp(predominant_faulting_style(n), 'ss') == 1
%                     if strcmp(invert.solve_for_correlation_length, 'yes') == 1
%                         a_as(n) = a_as_ss_mean + erfinv(m_trial(m_identifyer==10))*sqrt(2)*a_as_ss_std;   % Have to use error function trick to draw from a normal distribution, since we can't do a random walk with randn function
%                         a_dd(n) = a_dd_ss_mean + erfinv(m_trial(m_identifyer==11))*sqrt(2)*a_dd_ss_std;
%                     else
%                         a_as(n) =  1860 + 0.34*(slippingpatch_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
%                         a_dd(n) =  -390 + 0.44 * (slippingpatch_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002)
%                         
%                     end
%                 elseif strcmp(predominant_faulting_style(n), 'ds') == 1
%                     if strcmp(invert.solve_for_correlation_length, 'yes') == 1
%                         a_as(n) = a_as_ds_mean + erfinv(m_trial(m_identifyer==10))*sqrt(2)*a_as_ds_std;
%                         a_dd(n) = a_dd_ds_mean + erfinv(m_trial(m_identifyer==11))*sqrt(2)*a_dd_ds_std;
%                     else
%                         a_as(n) =  1100 + 0.31*(slippingpatch_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
%                         a_dd(n) =  580 + 0.35* (slippingpatch_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002)
%                     end
%                 end
%             end
%             
%    
%             % Recalculate scaled distances
%             [r_over_a, ~] = calc_scaled_dist( n_fault_strands_for_smoothing, disloc_model, a_as, a_dd, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master, along_strike_sep_dist, n_along_strike_patches_for_smoothing, n_down_dip_patches_for_smoothing,fault_strand_togetherness);
%             
%             % Recalculate sigma_s
%             sigma_s = calc_sigma_s( r_over_a(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)), H(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)));       % NOTE that this is actually a SEPARATE SIGMA S MATRIX FOR EACH FAULT STRAND, just stored in one big matrix
%             if strcmp(invert.add_correlation_matrix_stabiliser, 'yes') == 1
%                 sigma_s = sigma_s + diag( 0.01*ones(total_n_slip_patches,1));
%             end
%             sigma_s(:, onoffidentifyer==0) = [];            % slice columns out of master
%             sigma_s(onoffidentifyer==0, :) = [];            % slice rows out of master
%             
%             det_sigma_s = det(sigma_s);
%             inv_sigma_s = inv(sigma_s);     %Uncomment below for multiple fault strands
%             %det_sigma_s_keep(:,faultsizecount) = det_sigma_s;
%             
%                         % slip              % alpha                             % rake              %dip                                    % offset                                        % beta                   % circharms.
%             m_on = [onoffidentifyer; ones(n_fault_strands_for_smoothing,1); onoffidentifyer; ones(size(dip_initial,1),1); ones(size(offset_modelparameter_initial,1),1); ones(size(beta_initial,1),1); ones(sum(circharmparameters), 1)];
%             m_identifyer = m_identifyer_master;
%             m_identifyer(m_on==0) = [];  % Remove appropriate slip rows and rake rows from m_identifyer
%             
%             % Remove 'off' patches, and corresponding rakes, and G_temp
%             m_curr = m_initial;
%             %mfreeze = m_curr(m_on==0); % Keep the values of m_curr that are now off, for when these patches turn back on
%             m_trial = m_initial;        % this is necessary so all the matrices are the right lengths, plus circharm parameters are correct (since they don't change on every iteration)
%             m_initial(m_on==0) = [];
%             G_curr(:, onoffidentifyer==0) = [];
%             
%             % Clear not useful things
%             clear ch_coeffs_trial
%             clear ch_phi_trial
%             clear ch_center_trial
%             clear ch_coeffs_curr
%             clear ch_phi_curr
%             clear ch_center_curr
%             clear ch_coeffs_keep
%             clear ch_phi_keep
%             clear ch_center_keep
%             clear logprob_trial
%             clear logprob_curr
%             clear logprob_keep
%             clear ch_sigma_d
%             clear inv_ch_sigma_d
%             clear logL_curr_unscaled
%             clear logL_trial_unscaled
%             clear m_max_simulatedannealing
%             clear m_min_simulatedannealing
%         end
         
        % Show result
        figure; faults = disloc_model(:,onoffidentifyer==1); faults(6,:) = m_curr(m_on==1&m_identifyer_master==1)'; doplot3d(faults', 'jet'); colorbar
        title(['alpha^2 = ', num2str(alpha2_initial)]);
        disp('   ');
        disp('This is the simulated annealing result, which will be used for starting values.')
        disp('Click on the figure to continue...')
        disp('   ');
        %waitforbuttonpress;
        
        % Recalculate initial probabilities
        for ii = 1:nwalkers
            if strcmp(invert.smoothing, 'VK') == 1
                [logprior_curr(:,ii), ~] = calc_logprior_VK(m_initial(m_identifyer==1,ii), inv_sigma_s, alpha2_initial(:,ii), n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing);        % NOTE: simple means we don't calculate anything other than the exponent, because for non varying H and a, this is the same every time.
            end
            logL_curr(1,ii) = calc_loglikely(m_initial(m_identifyer==1,ii), d, G_curr, inv_sigma_d, offset_curr, m_initial(m_identifyer==6,ii));
            logposterior_curr(1,ii) =  log(sum( (2*pi*alpha2_initial(:,ii)).^(-n_slip_patches_on_each_fault_strand_for_smoothing/2))) + sum(logprior_curr(1,ii)) + logL_curr(:,ii);                      % if we hadn't taken logs, we'd multiply them. but we calculated the log probability, so we add them. log(AB) = log(A) + log(B).
        end
        if strcmp(invert.regularise_moment, 'yes') == 1
            M0_curr = sum((elastic_params.mu_okada) * spatial_model2column .* spatial_model3column .* (slip_initial(:,1)));
            M0_likelihood_curr = normpdf( M0_curr, data.seismic_moment, data.moment_std);
            logposterior_curr = logposterior_curr + log(M0_likelihood_curr);        % VK logposterior_curr
        elseif strcmp(invert.regularise_moment, 'no') == 1
            M0_likelihood_curr = [];
        end
        if strcmp(invert.solve_for_beta, 'yes') == 1
            logposterior_curr = logposterior_curr * (length(d)/2)*log(2*pi*m_initial(m_identifyer==6,ii)^2);
        end
        
        m_on_keep(:, 1) = m_on;
        m_keep = zeros(n_model_parameters, iterations, nwalkers);
        step_sizes(m_identifyer_master==1) = invert.step_size;
        step_sizes(m_identifyer_master==5) = abs(m_initial(m_identifyer==5)/10);