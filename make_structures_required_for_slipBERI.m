% This script makes the default structures for the slipBERI code
%
% For details on what each parameter is, please see the help for slipBERI

clear all

% first structure: 'fault'
fault = struct('fault_descriptor_file', 'name_of_text_file_with_fault_details_in.txt',...
                'fault_coordinate_unit', 'long/lat');       % or 'utm'

% second structure: 'data'
data = struct('InSAR_datafile', {{'name_of_InSAR_datafile.txt'}}, ...       % Name of InSAR text file, or 'none'
              'GPS_datafile_2d', 'name_of_GPS_2d_file.txt',...              % Name of GPS 2D text file, or 'none'
              'GPS_datafile_3d', 'name_of_GPS_3d_file.txt',...              % Name of GPS 3D text file, or 'none'
              'atolls_datafile', 'none',...                                 % atolls (or other dataset with only 'up' component)
              'InSAR_coordinate_unit', 'utm',...                            % 'utm' or 'long/lat'
              'GPS_coordinate_unit', 'long/lat',...                         % 'utm' or 'long/lat'
              'UTMzone', 10,...
              'origin', [133.815 35.368],...
              'EQ_epicenter', [35.368, 133.815, 5.7], ...     
              'seismic_moment', 2.818e+18, ...                              % Nm, e.g. from USGS. Only used if using moment regularisation
              'moment_std', 2.1310e+17,...                                  % Nm, e.g. from USGS. Only used if using moment regularisation
              'mindist', 0,...                                              % NOT FUNCTIONAL YET. choose whether you want to remove faults a certain distance from the fault or not
              'weight_InSAR', 1,...
              'weight_GPS', 1,...
              'weight_atolls', 1,...
              'varcovar_deets', {{'text_file_with_sill_nugget_range_for_covariance_calculation.txt'}},...
              'quadtree_n_points', {{'file_with_details_of_how_many_points_are_averaged_in_each_pixel.txt'}});       % if not downsampled using quadtree, this is just a text file of 1s, same length as InSAR data
          
% third structure: 'testing'
testing = struct('testing_mode', 'no',...
                 'making_model', 'name_of_saved_file_used_to_setup_test.mat',...
                 'var', 0,...
                 'add_noise', 0)  ;        
          
% fourth structure: 'invert'
invert = struct('inversion_type', 'bayesian', ...               % or 'least_squares'
                'smoothing', 'VK',...                           % or 'laplacian' or 'VK' or 'none'
                'smooth_across_fault_strands', 'yes',...        % or 'no'
                'iterations', 1000000, ...                      % only relevant for bayesian inversions
                'regularise_moment', 'no', ...                  % or 'yes'
                'solve_for_dip', 'no',...                       % NOT FULLY FUNCTIONAL
                'solve_for_correlation_length', 'no',...        % NOT FULLY FUNCTIONAL or 'yes'
                'solve_for_InSAR_offset', 'no',...              % or 'yes'
                'solve_for_fault_size', 'no', ...               % or 'yes'
                'slip_initial', 0.7, ...                        % meters
                'step_size', 0.05, ...                          % metres
                'alpha_initial', 0.01,...                
                'probability_target_initial', 0.001, ...
                'alpha_step_size', 0.01,...    
                'pad_edges_with_zeros', 'no',...                % NOT RECOMMENDED. or 'yes'        
                'variable_rake', 'yes',...                      % or 'no'
                'load_old_MCMC_chain', 'no',...                 % or name of old chain to be reloaded and inversion continued
                'solve_for_beta', 'no',...                      % or 'yes'
                'beta_initial', 1,...
                'beta_step_size', 0.01,...
                'simulated_annealing_start', 'no',...           % ONLY WORKS FOR VON KARMAN REGULARISATION. or 'yes'
                'ensemble_sampling', 'yes',...                  % NOT FULLY FUNCTIONAL
                'ensemble_move_style', 'stretch',...            % NOT FULLY FUNCTIONAL
                'ensemble_start', 'tight',...                   % NOT FULLY FUNCTIONAL
                'select_patches', 'no',...
                'circular_harmonics', 'no');       

% fifth structure: details on the priors
priors = struct('min_dip', 90,...                       % degrees
        'max_dip', 90,...                               % degrees
        'slip_prior', 'boxcar', ...            
        'min_slip', 0, ...                              % metres
        'max_slip', 50, ...                             % metres
        'min_rake', 150,...                             % degrees.
        'max_rake', 210,...                             % degrees.
        'alpha_prior', 'logarithmic',...
        'min_alpha2', 0.001,...
        'max_alpha2', 2,...
        'alphaflag', [],...
        'min_offset', -1,...                            % metres
        'max_offset', 1,...                             % metres
        'predominant_faulting_style', {{'ss'}},...      % 'ss' = strike-slip and 'ds' = dip-slip. These are used to calculate the von Karman autocorrelation lengths, using the appropriate equation for strike-slip or dip-slip
        'max_beta', 1,... 
        'min_beta', 1,...
        'min_circharm_coeffs', [], ...
        'max_circharm_coeffs', [], ...
        'min_circharm_phi', [], ...
        'max_circharm_phi', [], ...
        'min_circharm_center', [],...
        'max_circharm_center', []);                
        
% sixth structure: 'elastic_params'
elastic_params = struct('lambda', 3.23e10,...
                        'mu_okada', 3.23e10);

% seventh structure: 'disp', within which 'figures' is a nested structure of which figures you'd like
display = struct('plot_resolution_matrix', 'no',...     
                 'plotmean', 'no',...                        % other option: 'yes'
                 'plotmode', 'no',...                        % other option: 'yes'
                 'plotmedian', 'no',...                      % other option: 'no'              
                 'plotmaxlikely', 'yes',...                  % other option: 'yes'  % IN FUTURE VERSIONS THIS WILL BE DELETED, coz I've replaced it with MAP
                 'plotallslips', 'yes',...                   % other option: 'no'
                 'plotprob', 'yes',...                       % other option: 'no'
                 'plothists', 'plothistsall',...             % other option: 'plothistsample' (i.e. choose some at random to plot)
                 'plotsurfacedisp', 'yes',...                % other option: 'no'
                 'plotmarginalPDFs', 'no',...                % other option: 'no'
                 'plotMAP', 'yes',...                        % other option: 'no'
                 'plot3d', 'yes') ;                          % other option: 'no'
             
                
 % final
 housekeeping = struct('save_name', 'slip');      % your output file will be called '<savename>_parameters_used.mat'


% display all the structures so people actually look and check
disp('elastic_params'), disp(elastic_params)
disp('housekeeping ='); disp(housekeeping)
disp('testing ='), disp(testing)
disp('display ='), disp(display)
disp('fault ='); disp(fault)
disp('data ='); disp(data)
disp('invert ='); disp(invert)
disp('priors ='); disp(priors)

