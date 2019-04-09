% This script makes the default structures for the slipBERI code
%
% For details on what each parameter is, please see the help for slipBERI:
% https://docs.google.com/document/d/1cUXLRxN-oB8Q8kGOueq2c-Zxr3W1vDgWpGUw1MAEx5s/edit#


clear all

% 'fault' contains details about the fault model upon which slip will be inverted
fault = struct('fault_descriptor_file', 'name_of_text_file_with_fault_details_in.txt',...       % Name of fault text file
                'fault_coordinate_unit', 'long/lat');                                           % 'long/lat' or 'utm'

% 'data' contains details about the data which will be used in the inversion
data = struct('InSAR_datafile', {{'name_of_InSAR_datafile.txt'}}, ...       % {{'name_of_InSAR_text_file1.txt', 'name_of_InSAR_text_file2.txt'}}, or 'none'
              'InSAR_coordinate_unit', 'utm',...                            % 'long/lat' or 'utm'
              'varcovar_details', {{'text_file_with_sill_nugget_range_for_covariance_calculation.txt'}},...
              'quadtree_n_points', {{'text_file_with_details_of_how_many_points_are_averaged_in_each_pixel.txt'}},...       % If not downsampled using quadtree, this is just a text file of 1s, same length as InSAR data
              'GPS_datafile_2d', 'name_of_GPS_2d_file.txt',...              % Name of GPS 2D text file, or 'none'
              'GPS_datafile_3d', 'name_of_GPS_3d_file.txt',...              % Name of GPS 3D text file, or 'none'
              'GPS_coordinate_unit', 'long/lat',...                         % 'long/lat' or 'utm'
              'atolls_datafile', 'none',...                                 % Name of atolls text file, or 'none'
              'atolls_coordinate_unit', 'long/lat',...                      % atolls (or other dataset with only 'up' component)
              'weight_InSAR', 1,...                                         % Weighting of InSAR compared to GPS/atolls. If no weighting then this must be 1.
              'weight_GPS', 1,...                                           % Weighting of GPS compared to InSAR/atolls. If no weighting then this must be 1.
              'weight_atolls', 1,...                                        % Weighting of atolls compared to InSAR/GPS. If no weighting then this must be 1.
              'UTMzone', 10,...                                             % UTM zone of the earthquake. Ignore letter. Negative if in southern hemisphere.
              'origin', [133.815 35.368],...                                % A long/lat in the center of your area.
              'EQ_epicenter', [35.368, 133.815, 5.7], ...                   % Epicenter of earthquake.
              'seismic_moment', 2.818e+18, ...                              % Nm, e.g. from USGS. Only used if using moment regularisation
              'moment_std', 2.1310e+17);                                    % Nm, e.g. from USGS. Only used if using moment regularisation);                                          
                
          
% 'invert' contains details of how the inversion will be performed   
invert = struct('quickcheck', 'yes',...                         % 'yes' or 'no'
                'inversion_type', 'bayesian', ...               % 'bayesian' or 'least_squares'
                'iterations', 100000, ...                       % only relevant for bayesian inversions
                'smoothing', 'VK',...                           % 'VK' or 'laplacian' or 'none'
                'smooth_across_fault_strands', 'yes',...        % 'yes' or 'no'
                'slip_initial', 0.7, ...                        % meters
                'step_size', 0.05, ...                          % metres
                'variable_rake', 'yes',...                      % 'yes' or 'no'
                'solve_for_InSAR_offset', 'no',...              % 'yes' or 'no'
                'solve_for_InSAR_ramp', 'no',...                % 'yes' or 'no'
                'regularise_moment', 'no', ...                  % 'yes' or 'no'
                'alpha2_initial', 0.01,...    
                'alpha2_step_size', 0.01,... 
                'probability_target_initial', 0.001, ...
                'solve_for_beta', 'no',...                      % 'yes' or 'no'
                'beta_initial', 1,...
                'beta_step_size', 0,...
                'simulated_annealing_start', 'no', ...          % 'yes' or 'no'. Only works for von Karman regularisation.
                'solve_for_fault_size', 'no', ...               % 'yes' or 'no'
                'add_correlation_matrix_stabiliser', 'no',...     % use 'no' for standard von karman inversion, and 'yes' for trans-dimensional (i.e. when invert.solve_for_fault_size = 'yes')                
                'load_old_MCMC_chain', 'no');                 % 'no' or name of old chain to be reloaded and inversion continued
                        


% priors contains details about the prior information which will be used in the inversion
priors = struct('slip_prior', 'boxcar', ...                     % 'boxcar', 'gaussian' or 'logarithmic'
                'min_slip', 0, ...                              % metres
                'max_slip', 50, ...                             % metres
                'predominant_faulting_style', {{'ss'}},...      % 'ss' = strike-slip and 'ds' = dip-slip. These are used to calculate the von Karman autocorrelation lengths, using the appropriate equation for strike-slip or dip-slip
                'min_rake', 150,...                             % degrees.
                'max_rake', 210,...                             % degrees
                'alpha2_prior', 'logarithmic',...                % 'boxcar' or 'logarithmic'
                'min_alpha2', 0.001,...
                'max_alpha2', 2,...
                'alpha2flag', [],...                             % 'bothsame' or []
                'min_offset', -1,...                            % metres
                'max_offset', 1,...                             % metres
                'max_beta', 1,...
                'min_beta', 1,...
                'min_circharm_coeffs', [], ...
                'max_circharm_coeffs', [], ...
                'min_circharm_phi', [], ...
                'max_circharm_phi', [], ...
                'min_circharm_center', [],...
                'max_circharm_center', []);
        
% elastic_params contains values for Lame’s first and second parameters
elastic_params = struct('lambda', 3.23e10,...
                        'mu_okada', 3.23e10);

% display contains details of which parameters to plot at the end
display = struct('plot_resolution_matrix', 'no',...     
                 'plotmean', 'no',...                        % 'yes' or 'no'
                 'plotmode', 'no',...                        % 'yes' or 'no'
                 'plotmedian', 'no',...                      % 'yes' or 'no'
                 'plotMAP', 'yes',...                        % 'yes' or 'no'
                 'plotallslips', 'yes',...                   % 'yes' or 'no'
                 'plotprob', 'yes',...                       % 'yes' or 'no'
                 'plothists', 'plothistsall',...             % 'plothistall' or 'plothistsample'
                 'plotsurfacedisp', 'yes',...                % 'yes' or 'no'
                 'plotmarginalPDFs', 'yes',...               % 'yes' or 'no'
                 'plot3d', 'yes',...                         % 'yes' or 'no'
                 'calc_confidence', 'no') ;                  % 'yes' or 'no'
                         
 % housekeeping is the name of the output file
 housekeeping = struct('save_name', 'slip');      % your output file will be called '<savename>_parameters_used.mat'


% display all the structures so people actually look and check
disp('elastic_params'), disp(elastic_params)
disp('housekeeping ='); disp(housekeeping)
disp('display ='), disp(display)
disp('fault ='); disp(fault)
disp('data ='); disp(data)
disp('invert ='); disp(invert)
disp('priors ='); disp(priors)

