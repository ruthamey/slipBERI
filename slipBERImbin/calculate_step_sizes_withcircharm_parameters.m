function [ new_step_sizes ] = calculate_step_sizes_withcircharm_parameters( m, total_n_slip_patches, invert, G_curr, d,inv_sigma_d, inv_sigma_s, step_sizes, probability_target, G_ss_curr, G_ds_curr,n_fault_strands, n_fault_strands_for_smoothing, n_slip_patches_on_each_fault_strand, first_patch_in_strand, last_patch_in_strand, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master, det_sigma_s, L, M0_likelihood_before, elastic_params, spatial_model2, spatial_model3, data, m_identifyer_master, G_ss_LUT, G_ds_LUT, dip_LUT, n_data, fault_strand_identifyer, min_dip, max_dip, onoffidentifyer, n_slip_patches_ON_on_each_fault_strand_for_smoothing, n_slip_patches_on_each_fault_strand_for_smoothing, d_identifyer, patchx, patchz, disloc_model, n_down_dip_patches_for_smoothing, n_along_strike_patches_for_smoothing, along_strike_sep_dist, n_along_strike_patches, n_down_dip_patches, fault_strand_togetherness, H, first_patch_in_strand_for_smoothing,last_patch_in_strand_for_smoothing, fault_length_for_smoothing, fault_width_for_smoothing, InSAR_identifyer, n_harmonics, predominant_faulting_style, priors, n_InSAR_scenes, m_on)
% calculate_step_sizes_with_depth works out what step sizes should be used
% for each patch.
%
% An ideal probability perturbation is given, e.g. 0.1
% Then, each model parameter in turn is perturbated by half the current step size.
% The different in posterior probability resulting from this perturbation is
% calculated. The step size is then adjusted to ensure that the
% perturbation for each patch is equal to the ideal probability
% perturbation.
%
% This ideal probability perturbation should be updated through the 
% iteration, prior to this function, based upon the rejection ratio.
% If the rejection ratio is too high, the ideal probability target is made 
% smaller. This then persuades the step sizes to be made smaller,
%
% The final output of this function is a (1 x total_n_slip_patches) size
% matrix, where each patch has a value for step size
%
% Inputs:
%          m = [slip; alpha; rake]

% PLEASE SEE calculate_step_sizes_untidy for all the various ways of taking
% the product / sum of the probabilities, when I was getting confused about
% taking the ratio of log probabilities for multiple strands... I think
% this is right, so I've deleted the other versions. But they're in
% calculate_step_sizes_untidy if you'd like to check them.

% rmja 27-may-2015  wrote it
% rmja  2-sep-2015  changed it to perturb all patches and calculate average, not just a random patch
% rmja 29-sep-2015  changed it from calculate_step_sizes_with_depth to calculate_step_size for each patch, aiming towards an ideal probability perturbation
% rmja  7-oct-2015  add nan check
% rmja  3-nov-2015  correct for multiple fault strands
% rmja  6-jan-2016  added laplacian smoothing
% rmja  8-jan-2016  correct for M0 regularisation
% rmja  7-feb-2017  correct for solving for patches on and off

% Housekeeping
total_n_model_parameters = length(m);
%P_ratio = zeros(total_n_model_parameters, 1);
probability_after_perturbation = zeros(total_n_model_parameters, 1);
new_step_sizes = zeros(total_n_model_parameters, 1);
diff = zeros(total_n_model_parameters, 1);
slip_initial = m(m_identifyer_master ==1);
if strcmp(invert.solve_for_InSAR_offset, 'no') == 1
    n_InSAR_scenes = 0;
end
count = 1;                                                  % rake of slip patch 1 affects G(:,1). Then count to change each column of G in turn
if any(m_on==0)
    %m_on_synthetic_test = [onoffidentifyer; ones(n_fault_strands_for_smoothing,1); onoffidentifyer; ones(n_InSAR_scenes,1); zeros((n_harmonics*2-1+2),1) ];  % Slip AND rake AND circharmparameters are off, here         (slip, alpha, rake, dip, beta, offset, circharm coeff, circharm phi, cirdharm center)
    m_on_synthetic_test = [onoffidentifyer; ones(n_fault_strands_for_smoothing,1); onoffidentifyer; ones(n_InSAR_scenes,1); ones((n_harmonics*2-1+2),1) ];  % Slip AND rake are off, here         (slip, alpha, rake, dip, beta, offset, circharm coeff, circharm phi, cirdharm center)
else
   %m_on_synthetic_test = [ones(length(m_on)-(n_harmonics*2-1+2),1); zeros((n_harmonics*2-1+2),1) ];  % never solve for circharm parameters
   m_on_synthetic_test = [ones(length(m_on)-(n_harmonics*2-1+2),1); ones((n_harmonics*2-1+2),1) ];  % trying solving for circharm parameters
end
m_identifyer = m_identifyer_master;
m_identifyer(m_on_synthetic_test==0) = [];      % Turn off rake patches for slip patches that are off
% m_identifyer(m_identifyer==7) = [];
% m_identifyer(m_identifyer==8) = [];
% m_identifyer(m_identifyer==9) = [];
m_master = m;
m(m_on_synthetic_test==0) = [];
step_sizes_juston = step_sizes;
step_sizes_juston(m_on_synthetic_test==0) = [];
half_step_size = 0.5 * step_sizes_juston; % We'll be taking 0.5 the maximum step size in each case
% if strcmp(invert.smoothing, 'VK') || strcmp(invert.smoothing, 'laplacian') == 1
%     alpha_initial = m(total_n_slip_patches+1 : (total_n_slip_patches+n_fault_strands));
%     rake_start_number = n_fault_strands+1;                        % this is used to tell you where in your m matrix the rake comes in. m = [slip; alphas; rake], so we want to start rake at number (slip + number of alphas) and the number of alphas will be equal to the number of fault strands
% elseif strcmp(invert.smoothing, 'none') == 1;
%     rake_start_number = 1;                        % this is used to tell you where in your m matrix the rake comes in. m = [slip; alphas; rake], so we want to start rake at number (slip + number of alphas) and the number of alphas will be equal to the number of fault strand
% end
if strcmp(invert.select_patches, 'yes') ==1
    invert.solve_for_fault_size = 'no';
end

if strcmp(invert.smoothing, 'VK') || strcmp(invert.smoothing, 'laplacian')|| strcmp(invert.smoothing, 'minimumnorm') == 1
    [I, ~] = find( m_identifyer==2, 1);
    alpha_start_number = I;
    alpha_modelparameter_initial = m(m_identifyer== 2);
    if strcmp(priors.alpha_prior, 'logarithmic') == 1
        alpha_initial = 10.^alpha_modelparameter_initial;
    else
        alpha_initial = alpha_modelparameter_initial;
    end
    perturbed_alpha = alpha_initial;
end

if strcmp(invert.variable_rake, 'yes') == 1
    [I, ~] = find( m_identifyer==3, 1);
    rake_start_number = I;
    rake_number = rake_start_number;
    rake_initial = m( m_identifyer ==3);
end

if strcmp(invert.solve_for_dip, 'yes') ==1
    [I, ~] = find( m_identifyer==4, 1);
    dip_start_number = I;
    dip_number = dip_start_number;
else
    dip_start_number= [];
    dip_number= [];
end

if strcmp(invert.solve_for_InSAR_offset, 'yes') ==1
    [I, ~] = find( m_identifyer==5, 1);
    offset_start_number = I;
    offset_initialshort = m( m_identifyer == 5);
    offset_initial = zeros(length(d),1);                        % make the offset vector the same length as the d vector...
    for pq = 1: length(offset_initialshort)
        offset_initial(InSAR_identifyer==pq) = offset_initialshort(pq);        %... but only add offset to InSAR. And add the one offset value to ALL the InSAR
    end
else
    offset_initial = zeros(length(d),1);
    offset_start_number = [];
end

    beta_initial = 1;
    perturbed_beta = 1;

    [I, ~] = find( m_identifyer_master==6, 1);
    beta_start_number = I;
    if strcmp(invert.solve_for_beta, 'yes') == 1           % if we're solving for beta
        beta_initial = m(m_identifyer_master == 6);
        perturbed_beta = beta_initial;
    end
    
    if strcmp(invert.solve_for_fault_size, 'yes') == 1
       G_curr_master = G_curr;
       slip_initial(m_on==0) = [];
       G_curr(:, onoffidentifyer==0) = [];
       G_ss_curr(:, onoffidentifyer==0) = [];
       G_ss_LUT(:, onoffidentifyer==0) = [];
       G_ds_curr(:, onoffidentifyer==0) = [];
       G_ds_LUT(:, onoffidentifyer==0) = [];
       inv_sigma_s(:, onoffidentifyer==0) = [];
       inv_sigma_s(onoffidentifyer==0, :) = []; 
    end

% Calculate probability of the initial given slip distribution
if strcmp(invert.inversion_type, 'bayesian') == 1 && strcmp(invert.smoothing, 'none') == 1
            %d_hat_initial = G_curr * slip_initial;
            probability_before_perturbation = calc_loglikely(slip_initial, d, G_curr, inv_sigma_d, offset_initial, beta_initial);        %  not exp((-1/2) * ( (d - d_hat).' * inv_sigma_d * (d - d_hat) ));
            probability_before_perturbation = (-length(d)/2)*log((2*pi*perturbed_beta^2)) + probability_before_perturbation;
elseif strcmp(invert.smoothing, 'VK') == 1
            exponent_VK_before = calc_logprior_VK( slip_initial, inv_sigma_s, alpha_initial, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing);
            logL_before = calc_loglikely(slip_initial, d, G_curr, inv_sigma_d, offset_initial, beta_initial);
            probability_before_perturbation = sum((-n_slip_patches_ON_on_each_fault_strand_for_smoothing/2).*log(2*pi*alpha_initial.^2)) + (-length(d)/2)*log((2*pi*beta_initial^2) ) + sum(exponent_VK_before) + logL_before;          % matrix matrix. need to keep as a matrix coz the constants out the front are different for multiple faults strands. dont' need to calculate them as long as we work them out separately
            % this is actually the LOG probability. sum over all slip patches.
elseif strcmp(invert.smoothing, 'laplacian') == 1
            exponent_laplacian_before = calc_logprior_laplacian( slip_initial, L, alpha_initial, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing, onoffidentifyer);
            logL_before = calc_loglikely(slip_initial, d, G_curr, inv_sigma_d, offset_initial, beta_initial);
            probability_before_perturbation = sum((-n_slip_patches_ON_on_each_fault_strand_for_smoothing/2).*log(2*pi*alpha_initial.^2)) + (-length(d)/2)*log((2*pi*beta_initial^2)) + sum(exponent_laplacian_before) + logL_before;
elseif strcmp(invert.smoothing, 'minimumnorm') == 1
            %logprior_before =  alpha_initial * (slip_initial' * slip_initial);
            for q = 1:n_fault_strands_for_smoothing
                s_curr = slip_initial(first_patch_in_strand_for_smoothing(q):last_patch_in_strand_for_smoothing(q));
                logprior_before(q,1) =  ( -1/(2*alpha_initial(q).^2) * sum((s_curr'*s_curr)));
            end
            logL_before = calc_loglikely(slip_initial, d, G_curr, inv_sigma_d, offset_initial, beta_initial);    
            probability_before_perturbation = -(length(d)/2)*log((2*pi*beta_initial^2)) + sum(log(2*pi*alpha_initial)) + logL_before + sum(logprior_before);          % matrix matrix. need to keep as a matrix coz the constants out the front are different for multiple faults strands. dont' need to calculate them as long as we work them out separately
end

if strcmp( invert.regularise_moment, 'yes') == 1
   probability_before_perturbation =  probability_before_perturbation + log(M0_likelihood_before);
end


%allthei = 1: total_n_model_parameters;
number_of_on_model_parameters = length(m_identifyer);

tofindtherighti = cumsum(m_on_synthetic_test);

for i = 1: number_of_on_model_parameters
    
     % Perturb model parameters in turn (each slip patch, then alpha...)  
     perturbed_m = m;
     perturbed_G = G_curr;
     perturbed_m(i) = m(i) + half_step_size(i); 
     mi = find(tofindtherighti==i,1, 'first');       % we're looping through the on model parameters. this tells us it correct position in all the model parameters, so we can update correct place in diff matrix.
     perturbed_slip = perturbed_m(m_identifyer == 1);
     perturbed_offset = offset_initial;

      if strcmp(invert.smoothing, 'VK') || strcmp(invert.smoothing, 'laplacian') || strcmp(invert.smoothing, 'minimumnorm') == 1
         if strcmp(priors.alpha_prior, 'logarithmic') == 1
                perturbed_alpha = 10.^perturbed_m(m_identifyer == 2);
         else
             perturbed_alpha = perturbed_m( m_identifyer == 2);          % either update it, or set it back to the original value, coz we do dip/rake afterwards. also have to give all of the alphas, even tho we've just perturbed one of them, or maybe not even any of them
         end
      end
     
      if m_identifyer(i) == 1       %solving for slip

      elseif m_identifyer(i) == 2   %solving for           already done this, coz you need to either reset alpha to its unperturbed value / update it, coz alpha is used in EVERY calculation.
 
            %disp('hi!')
            %disp(num2str(i));
            %keyboard
            % perturbed_alpha = perturbed_m( m_identifyer == 2);  %ALREADY DONE!    % all of the alphas, if there are multiple slip strands
            if strcmp(priors.alphaflag, 'bothsame') == 1
               perturbed_alphamodelparameter = perturbed_m(i)*ones(n_fault_strands_for_smoothing,1);
               if strcmp(priors.alpha_prior, 'logarithmic') == 1
                    perturbed_alpha = 10.^perturbed_alphamodelparameter;
               else
                    perturbed_alpha = perturbed_alphamodelparameter;          % either update it, or set it back to the original value, coz we do dip/rake afterwards. also have to give all of the alphas, even tho we've just perturbed one of them, or maybe not even any of them
               end
            end
            
      elseif m_identifyer(i) == 3   %solving for rake
            perturbed_rake = perturbed_m( rake_number);  % select the rake corresponding
            %theta = 180 - perturbed_rake;
            theta = deg2rad(perturbed_rake);
            perturbed_G(:,count) = cos(theta) * G_ss_curr(:,count) + sin(theta) * G_ds_curr(:,count);    % just need to update the appropriate row. can use same G since we've not changed the fault
            count = count + 1;
            rake_number = rake_number +1;    

      elseif m_identifyer(i)  == 4  %solving for dip
        
          %G_ss_perturbed = G_ss_curr;       % we just need to update G_ss and G_ds for the columns corresponding to the patches on the fault strand that we've changed the dip of. if you get me.
          %G_ds_perturbed = G_ds_curr;       % as in, if we only perturb dip on fault 1, we want to change the G_ss values relating to fault 1, but we can leave the rest the same.
          
          if i == dip_start_number   % start counter again
                count = 1;
          end
          perturbed_dip = perturbed_m( dip_number );        % just select the dip that's been perturbed - just need to update G relating to that fault strand
          
          if perturbed_dip > max_dip
              perturbed_dip = max_dip;
          elseif perturbed_dip < min_dip
              perturbed_dip = min_dip;
          end
          
          %theta = 180 - rake_initial;
          theta = degtorad(rake_initial);
              
           % JUST updated G_ss and G_ds for the appropriate columns, relating to the fault that we've changed dip on. that's why I'm using count.
           G_ss_temp(:, first_patch_in_strand(count):last_patch_in_strand(count)) = interp3(1:n_slip_patches_on_each_fault_strand(count), 1:n_data, dip_LUT, G_ss_LUT(:,first_patch_in_strand(count):last_patch_in_strand(count),:), 1:n_slip_patches_on_each_fault_strand(count), 1:n_data, perturbed_dip);
           G_ds_temp(:, first_patch_in_strand(count):last_patch_in_strand(count)) = interp3(1:n_slip_patches_on_each_fault_strand(count), 1:n_data, dip_LUT, G_ds_LUT(:,first_patch_in_strand(count):last_patch_in_strand(count),:), 1:n_slip_patches_on_each_fault_strand(count), 1:n_data, perturbed_dip);

            
            %perturbed_G(:,count) = cos(theta) * G_ss_temp(:,count) + sin(theta) * G_ds_temp(:,count);    % just need to update the appropriate row
                  
           for r = first_patch_in_strand(count) : last_patch_in_strand(count)      % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake. only for slip patches on the fault that's had its dip changed.
                perturbed_G(:,r) = cos(theta(r)) * G_ss_temp(:,r) + sin(theta(r)) * G_ds_temp(:,r);
           end 
            
           count = count + 1;
           rake_number = rake_number +1;
           dip_number = dip_number + 1;
           
      elseif m_identifyer(i) == 5
          
          perturbed_offset = zeros(length(d),1);
          if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
              perturbed_offsetshort = perturbed_m( m_identifyer == 5);
              for pq = 1:length(perturbed_offsetshort)
                  perturbed_offset(InSAR_identifyer==pq) = perturbed_offsetshort(pq);
              end
          end          
          count = count + 1;      
      
      elseif m_identifyer(i) == 6
          
          perturbed_beta = perturbed_m(i);
         
      
      elseif m_identifyer(i) == 7 || m_identifyer(i) == 8 || m_identifyer(i) == 9     % if we're changing the size of the slipping area, need to recalculate which patches are on
         
              circharm_coeffs_trial = perturbed_m(m_identifyer == 7);
              circharm_phi_trial = perturbed_m(m_identifyer == 8);
              circharm_center_trial = perturbed_m(m_identifyer == 9);

              for k = 1 : n_harmonics
                  if circharm_coeffs_trial(k) < 0       % Make sure magnitude of circular harmonics stays positive. stay positive, you hear me?!
                      circharm_coeffs_trial(k) =  abs(circharm_coeffs_trial(k));
                  end
              end

              if circharm_center_trial(1) < 0       % Make sure center of slipping patch stays within fault
                  circharm_center_trial(1) =  abs(circharm_center_trial(1));
              elseif circharm_center_trial(1) > fault_length_for_smoothing
                  circharm_center_trial(1) = 2*fault_length_for_smoothing - circharm_center_trial(1);
              end

              if circharm_center_trial(2) < 0       % Make sure center of slipping patch stays within fault
                  circharm_center_trial(2) =  abs(circharm_center_trial(2));
              elseif circharm_center_trial(2) > fault_width_for_smoothing
                  circharm_center_trial(2) = 2*fault_width_for_smoothing - circharm_center_trial(2);
              end

              % calc circular harmonics
              [circx,circz]=circharm(circharm_coeffs_trial,circharm_phi_trial,0);     % set to '1' if want to plot

              % sort out onoffidentifyer
              newonoffidentifyer = inpolygon(patchx, patchz, (circx+circharm_center_trial(1)), (circz+circharm_center_trial(2)))';  % new onoffidentifyer is the size of all the slip patches
              %newpatchestoturnoff = newonoffidentifyer(onoffidentifyer==1);
              perturbed_slip = perturbed_m(m_identifyer_master==1);
              perturbed_slip(newonoffidentifyer==0)= [];
              perturbed_G = G_curr_master;
              perturbed_G(:, newonoffidentifyer==0) = [];
              
              
                        tippytop = min(disloc_model(8,newonoffidentifyer==1));     % this finds the minimum top depth of any patches that are on from the matrix disloc model 
                        verybottom = max(disloc_model(8,newonoffidentifyer==1));   % this finds the maximum top depth of any patches that are on from the matrix disloc model
                        slippingpatch_width_for_smoothing =  verybottom-tippytop;

                        farleft = min(disloc_model(1,newonoffidentifyer==1));
                        farright = max(disloc_model(1,newonoffidentifyer==1));
                        slippingpatch_length_for_smoothing = abs(farleft-farright);           % abs in case in a synthetic test you've defined +x and -x

                        if length(slippingpatch_width_for_smoothing) == 0
                           slippingpatch_width_for_smoothing = 0; 
                        end

                        if length(slippingpatch_length_for_smoothing) == 0
                            slippingpatch_length_for_smoothing = 0;
                        end

                        % work out how many patches in each strand are 'on'
                        for n = 1:n_fault_strands_for_smoothing
                            n_slip_patches_ON_on_each_fault_strand_for_smoothing(n,1) = sum(newonoffidentifyer(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),1));
                        end
                            first_patch_in_strand_for_smoothing(1) = 0;
                            first_patch_in_strand_for_smoothing(2:(n_fault_strands_for_smoothing+1),1) = cumsum(n_slip_patches_ON_on_each_fault_strand_for_smoothing);
                            first_patch_in_strand_for_smoothing = first_patch_in_strand_for_smoothing+1;
                            first_patch_in_strand_for_smoothing(n_fault_strands_for_smoothing+1) = []; 
                            last_patch_in_strand_for_smoothing = cumsum(n_slip_patches_ON_on_each_fault_strand_for_smoothing); 


                        % recalculate sigma_s, det_sigma_s and inv_sigma_s if doing von Karman inver sion

                            if strcmp(invert.smoothing, 'VK') == 1
                                if length(slippingpatch_width_for_smoothing) ~= 1
                                    keyboard
                                end
                                for n = 1 : n_fault_strands_for_smoothing 
                                    if strcmp(predominant_faulting_style(n), 'ss') == 1
                                        a_as(n) =  1860 + 0.34*(slippingpatch_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
                                        a_dd(n) =  -390 + 0.44 * (slippingpatch_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002) 
                                    elseif strcmp(predominant_faulting_style(n), 'ds') == 1 
                                        a_as(n) =  1100 + 0.31*(slippingpatch_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
                                        a_dd(n) =  580 + 0.35* (slippingpatch_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002) 
                                    end
                                end
                                [r_over_a, ~] = calc_scaled_dist( n_fault_strands_for_smoothing, disloc_model, a_as, a_dd, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master, along_strike_sep_dist, n_along_strike_patches, n_down_dip_patches, fault_strand_togetherness);
                                sigma_s = [];                          
                                        for n = 1: n_fault_strands_for_smoothing
                                            sigma_s_temp = calc_sigma_s( r_over_a(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)), H(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)));       % NOTE that this is actually a SEPARATE SIGMA S MATRIX FOR EACH FAULT STRAND, just stored in one big matrix   
                                            sigma_s = blkdiag(sigma_s,sigma_s_temp);
                                            sigma_s_temp = [];
                                        end
                                sigma_s_master = sigma_s;       % sigma_s_master is updated because faultlength and faultwidth changed, so we had to recalculate everything
                                                                % sigma_s only contains the rows and columns that correspond to patches that are 'on'.
                                sigma_s(:, newonoffidentifyer==0) = [];            % slice columns out of master
                                sigma_s(newonoffidentifyer==0, :) = [];            % slice rows out of master

                                det_sigma_s = [];
                                inv_sigma_s = [];
                                    for n = 1: n_fault_strands_for_smoothing            % NOTE that this is actually a SEPARATE invSIGMA S MATRIX FOR EACH FAULT PATCH, just stored in one big matrix, in the third dimension for each fault strand
                                        det_sigma_s(n,1) = det(sigma_s(first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n),first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n)));       
                                        inv_sigma_s_temp(:, :) = inv(sigma_s(first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n),first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n)));
                                        inv_sigma_s = blkdiag(inv_sigma_s, inv_sigma_s_temp);
                                        inv_sigma_s_temp = [];
                                    end

                               %det_sigma_s_ratio = prod((det_sigma_s_keep(:,faultsizecount)/det_sigma_s_keep(:,faultsizecount-1))^(-0.5));   % But we only need to calculate this for the change between n patches being on and m patches being on. once a trial is accepted, m patches will remain on for the next 1000 iterations. so once a trial is accepted we need to set this back to 1      
                               m_identifyer = m_identifyer_master;
                               if any(m_on==0)
                                   m_on_synthetic_test = [newonoffidentifyer; ones(n_fault_strands_for_smoothing,1); newonoffidentifyer; ones(n_InSAR_scenes,1); ones((n_harmonics*2-1+2),1) ];  % Slip AND rake AND circharmparameters are off, here         (slip, alpha, rake, dip, beta, offset, circharm coeff, circharm phi, cirdharm center)
                               else
                                   m_on_synthetic_test = [ones(length(m_on)-(n_harmonics*2-1+2),1); ones((n_harmonics*2-1+2),1) ];  % never solve for circharm parameters
                               end
                               m_identifyer(m_on_synthetic_test==0) = [];
                               
                            end
      end
     
%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%% Calculate probabilities and ratios %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%
        if strcmp(invert.inversion_type, 'bayesian') == 1 && strcmp(invert.smoothing, 'none') == 1
                %d_hat = perturbed_G * perturbed_slip;
                probability_after_perturbation = calc_loglikely(perturbed_slip, d, perturbed_G, inv_sigma_d, perturbed_offset, perturbed_beta);        %  not exp((-1/2) * ( (d - d_hat).' * inv_sigma_d * (d - d_hat) ));
                probability_after_perturbation = (-length(d)/2)*log((2*pi*perturbed_beta^2)) + probability_after_perturbation;
        elseif strcmp(invert.smoothing, 'VK') == 1
                exponent_VK_after = calc_logprior_VK( perturbed_slip, inv_sigma_s, perturbed_alpha, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing);
                logL_after = calc_loglikely(perturbed_slip, d, perturbed_G, inv_sigma_d, perturbed_offset, perturbed_beta);
                probability_after_perturbation = sum((-n_slip_patches_ON_on_each_fault_strand_for_smoothing/2).*log(2*pi*perturbed_alpha.^2)) + (-length(d)/2)*log((2*pi*perturbed_beta^2)) + sum(exponent_VK_after) + logL_after;     % this is actually the LOG probability. sum over all slip patches.    
        elseif strcmp(invert.smoothing, 'laplacian') == 1
                exponent_laplacian_after = calc_logprior_laplacian( perturbed_slip, L, perturbed_alpha, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing);
                logL_after = calc_loglikely(perturbed_slip, d, perturbed_G, inv_sigma_d, perturbed_offset, perturbed_beta);
                probability_after_perturbation = sum((-n_slip_patches_ON_on_each_fault_strand_for_smoothing/2).*log(2*pi*perturbed_alpha.^2)) + (-length(d)/2)*log((2*pi*perturbed_beta^2)) + sum(exponent_laplacian_after) + logL_after;              
        elseif strcmp(invert.smoothing, 'minimumnorm') == 1
            logprior_after =  perturbed_alpha * (perturbed_slip' * perturbed_slip);
            for q = 1:n_fault_strands_for_smoothing
                s_curr = perturbed_slip(first_patch_in_strand_for_smoothing(q):last_patch_in_strand_for_smoothing(q));
                logprior_after(q,1) =  ( -1/(2*perturbed_alpha(q).^2) * sum((s_curr'*s_curr)));
            end
            logL_after = calc_loglikely(perturbed_slip, d, perturbed_G, inv_sigma_d, perturbed_offset, perturbed_beta);    
            probability_after_perturbation = (-length(d)/2)*log((2*pi*perturbed_beta^2)) + sum(log(2*pi*perturbed_alpha)) + logL_after + sum(logprior_after);          % matrix matrix. need to keep as a matrix coz the constants out the front are different for multiple faults strands. dont' need to calculate them as long as we work them out separately
        end

       if strcmp( invert.regularise_moment, 'yes') == 1
            M0_after = sum((elastic_params.mu_okada) * (reshape(spatial_model2, total_n_slip_patches,1) .* reshape(spatial_model3, total_n_slip_patches,1)) .* (perturbed_slip));
            M0_likelihood_after = normpdf( M0_after, data.seismic_moment, data.moment_std);
            probability_after_perturbation =  probability_after_perturbation + log(M0_likelihood_after);
        end
        

        
        % make sure P ratio is always positive 
        if abs(probability_after_perturbation) > abs(probability_before_perturbation)
            P_ratio = 1-(probability_before_perturbation / probability_after_perturbation);    % save a ratio for each patch
        else
            P_ratio = 1-(probability_after_perturbation / probability_before_perturbation);            % save a ratio for each patch
        end
        
     % Calculate how different the ratio is to the target
     diff(mi) = probability_target - P_ratio;
        
        
     % Change the step size, so that the ratio of probability before and after is equal to the ideal_probability_perturbation ( nicked this formula from Andy)
     if diff(mi) > 0                                                                 % if it's making too little difference....
         %new_step_sizes(i,1) = step_sizes(i) * exp( diff(mi)/ (1-probability_target)/2); % ... make the step size bigger. exp(positive_number).
         new_step_sizes(mi,1) = step_sizes(mi) * exp( diff(mi)/ (probability_target*16)); % ... make the step size bigger. exp(positive_number).
     elseif diff(mi) < 0                                                             % if it's making too much difference....
        if strcmp(invert.regularise_moment, 'yes') == 1
            new_step_sizes(mi,1) = step_sizes(mi) * exp( diff(mi)/ abs(invert.probability_target_initial-probability_target)/8); % ... make the step size smaller. exp(negative_number). it's just an empirical formula, you can play around with the number if you'd like it to get smaller faster.
        elseif strcmp(invert.smoothing, 'none') || strcmp(invert.smoothing, 'laplacian') == 1
            new_step_sizes(mi,1) = step_sizes(mi) * exp( diff(mi)/abs(invert.probability_target_initial-probability_target)/8);
        else
           new_step_sizes(mi,1) = step_sizes(mi) * exp( diff(mi)/abs(invert.probability_target_initial-probability_target)/8);  % I think this is the best one, but not for laplacian or unsmoothed, for some reason. 
        end
     elseif diff(mi) == 0
         new_step_sizes(mi,1) = step_sizes(mi);
     end

     
end

% % % Plot the P-ratio
% figure('position', [300, 300, 1500, 900]);                    % slip, alpha, rake, dip, offset, beta 
% bar(1:total_n_slip_patches, P_ratio(m_identifyer==1), 'r');
% hold on
% bar(alpha_start_number, P_ratio(m_identifyer==2), 'b');
% bar(rake_start_number:(rake_start_number+total_n_slip_patches-1),P_ratio(m_identifyer==3), 'm');
% bar(offset_start_number, P_ratio(m_identifyer==5), 'k');
% bar(beta_start_number, P_ratio(m_identifyer==6), 'g');
% legend('slip', 'alpha', 'rake', 'offset', 'beta', 'Location', 'northeastoutside');
% plot([1,length(m)], [probability_target,probability_target]);
% title('P ratio')
% text(ceil(total_n_slip_patches/2),(probability_target+(0.1*probability_target)),'P ratio too big - step sizes will be reduced')
% text(ceil(total_n_slip_patches/2),(probability_target-(0.1*probability_target)),'P ratio too small - step sizes will be increased')
% 
% % % % % % % % % % % % % % % % % Plot the diff
% figure('position', [300, 300, 1500, 900]);                    % slip, alpha, rake, dip, offset, beta
%  bar(1:sum(m_identifyer_master==1), diff(m_identifyer_master==1), 'r');
%   hold on
%    if strcmp(invert.smoothing, 'none') == 0
%       bar(find(m_identifyer_master==2,1, 'first'):find(m_identifyer_master==2,1, 'last'), diff(m_identifyer_master==2), 'b');
%    end
%  bar(find(m_identifyer_master==3,1, 'first'):find(m_identifyer_master==3,1, 'last'), diff(m_identifyer_master==3), 'm');
%  if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
%     bar(find(m_identifyer_master==5,1, 'first'):find(m_identifyer_master==5,1, 'last'), diff(m_identifyer_master==5), 'g');
%  end
%  if strcmp(invert.solve_for_beta, 'yes') == 1
%     bar(beta_start_number, diff(m_identifyer_master==6), 'k');
%  end
%  legend('slip', 'alpha', 'rake', 'offset', 'beta', 'Location', 'northeastoutside');
%  plot([1,length(m)], [probability_target,probability_target]);
%  title('diff')
%  text(ceil(total_n_slip_patches/2),(max(diff)+0.5*max(diff)),'diff positive = step sizes will be increased')
%  text(ceil(total_n_slip_patches/2),(0.1*min(diff)),'diff negative = step sizes will be reduced')

% %  % %Prevent step sizes from being too big or too small
%new_step_sizes(new_step_sizes< 1e-20) = 1e-20;    
%new_step_sizes(new_step_sizes> 5) = 5;

% Special line to make sure offset step size doesn't get too big - that can mess everything up
% if new_step_sizes(m_identifyer==5) > 0.5
%     new_step_sizes(m_identifyer==5) = 0.5;  % don't let offset step size get bigger than 0.5
% end

% Special line to make sure slip step size doesn't get too big
if new_step_sizes(m_identifyer_master==1) < 0.001
    new_step_sizes(m_identifyer_master==1) = 0.001;
end

new_step_sizes(m_on_synthetic_test==0) = step_sizes(m_on_synthetic_test==0);    %don't change step sizes for patches that are off. THIS IS DONE ABOVE SINCE diff(0) = old step sizes.

% Special line to make sure rake step size doesn't get too small
if new_step_sizes(m_identifyer_master==3) < 0.1
    new_step_sizes(m_identifyer_master==3) = 0.1; 
end

new_step_sizes = abs(new_step_sizes);   % always positive

%if strcmp(invert.select_patches, 'yes') ==1;
%    invert.solve_for_fault_size = 'yes';
%end

end

