%remove_burn_in(burn_in_remove_number)
% Function to remove burn in, as dictated by the user

m_keep(:,1:burn_in_remove_number/nwalkers,:)=[];
logposterior_keep(:,1:burn_in_remove_number) = [];
if strcmp(invert.regularise_moment, 'yes') == 1
    M0_keep(:,1:burn_in_remove_number) = [];
    M0_likelihood_keep(:,1:burn_in_remove_number) = [];
end
logL_keep(:,1:burn_in_remove_number)=[];


% Separate each parameter in m_keep. These are used in display_result
slip_keep = m_keep(m_identifyer_master==1,:,:);
if strcmp(priors.alpha2_prior, 'logarithmic') ==1
    alpha2_keep = 10.^(m_keep(m_identifyer_master==2,:,:));
else
    alpha2_keep = m_keep(m_identifyer_master==2,:,:);
end
rake_keep = m_keep(m_identifyer_master==3,:,:);
dip_keep = m_keep(m_identifyer_master==4,:,:);
offset_keep = m_keep(m_identifyer_master==5,:,:);
beta_keep = m_keep(m_identifyer_master==6,:,:);
circharm_coeffs_keep = m_keep(m_identifyer_master==7,:,:);
circharm_phi_keep = m_keep(m_identifyer_master==8,:,:);
circharm_center_keep = m_keep(m_identifyer_master==9,:,:);
ramp_keep = m_keep(m_identifyer_master==12,:,:);

disp('Resaving now burn in has been removed.')
savename = [housekeeping.save_name, '_', num2str(n_down_dip_patches_for_smoothing(1)), 'x', num2str(n_along_strike_patches_for_smoothing(1)), '_', invert.smoothing, 'smooth_', invert.solve_for_dip, 'dip_', num2str(invert.iterations), '_', invert.regularise_moment, 'M0reg_', priors.slip_prior, '_', invert.solve_for_fault_size, 'patchesonoff'];            
save(savename, '-v7.3');

disp('Displaying result again...')
disp('...')
disp('..')
disp('.')
disp('Keyboard mode now. To terminate keyboard mode and end the slipBERI function, type ''dbcont'' and press Enter')