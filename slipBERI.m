 function [  ] = slipBERI( fault, data, testing, invert, priors, elastic_params, display, housekeeping )
%
% slipBERI is a code to invert for earthquake slip. 
%
% The inversion options are:
% - Least-squares. This includes NO smoothing and tikhonov smoothing.
% - Bayesian, with no smoothing, or Laplacian smoothing, or von Karman smoothing.
%
% The von Karman prior option allows for incorporating the evidence from
% several sources that fault roughness, and slip, shows fractal properties.
% The von Karman distribution is used as the prior slip distribution, and 
% is incorporated  using with Bayes rule. Monte Carlo Markov chain sampling 
% is used to sample the distribution.
%
% The inputs are listed below. I recommend running 'make_default_structure.m' 
% and then changing the appropriate values.
% Structures you will definitely have to change: fault, data, invert
% Structures you may have to change: display, testing (if you're testing), housekeeping
% Structures you're unlikely to change: elastic_params
% Values are changed by identifying what you'd like to change (e.g.
% 'min_slip') and which structure it's in (e.g. 'invert') and then typing into
% the terminal: 
%                structure.name_of_variable = <new_value> 
% For example:
%                   invert.min_slip = 0
% NOTE that if the variable you are changing is a string, then obviously it
% must have quotation marks around it.
% For example:
%                housekeeping.savename = 'slipsolution'
%
% This code should be run in the same folder as your datafiles. Or the
% datafiles must be added to path before starting this code.
%
%
%
% *********** % ************ % ********** % ********** % ********** % *****
% The inputs are seven structures containing this information:
%
% INPUTS: 
%     'fault' is a structure relating to the description of the fault upon which you'd like to invert:
%           fault_descriptor_file = String. Name of text file of: [strike, dip, rake, up dip projection of CENTRE of fault plane, up dip projection of CENTRE of fault plane, length (km), top depth (km), bottom depth (km), patches along strike, patches down dip, smoothingidentifyer]. where fault strands with the same smoothingidentifyer will be smoothed across, if you select to smooth across fault strands at all. If you please.      
%           fault_coordinate_unit = String. 'utm' or 'long/lat' to indicate whether the coordinates used in fault.fault_desriptor_file are in given in UTM or long/lat.
% 
%     'data' is a structure containing details of the data that you're inverting:
%           InSAR_datafile = Text file of [x (long), y (lat), observed LOS displacement, (E, N, up) component of LOS vector]. Make sure you've downsampled and added pointing vector. LOS vector convention is positive looking at the satellite. If  you have no InSAR data this MUST be 'none'. Locations translated to UTMx or a local latitude, if spanning multiple UTM zones.
%           GPS_datafile_2d = GPS file of 2d surface displacement. Set upto be the same as InSAR. Text file of 7 columns [x (long), y (lat), observed LOS displacement (m), LOS vector E, LOS vecor N, LOS vector up, sigma] where sigma is the standard deviation (mm for GPS, m for atolls). Each GPS station will have two rows: one for E LOS vector, one for N LOS vector. The sigma must be for E, N or up appropriately depending on which value you have. If  you have no GPS data this MUST be 'none'. Locations translated to UTMx or a local latitude, if spanning multiple UTM zones.
%           GPS_datafile_3d = GPS file of 3d surface displacement. Set upto be the same as InSAR. Text file of 7 columns [x (long), y (lat), observed LOS displacement (m), LOS vector E, LOS vecor N, LOS vector up, sigma] where sigma is the standard deviation (mm for GPS, m for atolls). Each GPS station will have three rows: one for E LOS vector, one for N LOS vector, one for up LOS vector. The sigma must be for E, N or up appropriately depending on which value you have. If  you have no GPS data this MUST be 'none'. Locations translated to UTMx or a local latitude, if spanning multiple UTM zones.
%           GPS_coordinate_unit = String. 'lat/long' or 'UTM'. Everything will be converted to UTM as necessary.
%           EQ_epicenter = Number. Matrix of [lat, long, depth]. Depth in km.
%           seismic_moment = Number. Only used if invert.regularise_moment == 'yes'. Estimate of seismic moment e.g. from USGS in Nm
%           moment_std = Number. Only used if invert.regularise_moment == 'yes'. Estimate of seismic moment e.g. variation of USGS solutions
%           mindist = Number. The distance from the fault from which data points will be removed - to minimise errors associated from incorrect fault geometry. NOT WORKING YET.
%           weight_InSAR = Number. Relative weighting of InSAR data. If not weighting datasets then this must be 1. If this is a number other than 1 then the InSAR var-covar matrix will be divided by this value.
%           weight_GPS = Number. Relative weighting of GPS data. If not weighting datasets then this must be 1. If this is a number other than 1 then the GPS var-covar matrix will be divided by this value.
%           weight_atolls = Number. Relative weighting of atoll dataset.
%           varcovar_deets = String. this is a text file with the [sill, variogram_range] from the covariogram from an undeforming region. Can be calculated in advance using andy's code variogram.m; In older versions this was: Name of a text file containing the outputsfrom cvdcalc: [maxium covariance in metres^2, alpha from expcos function, beta from expcos function, eb_r from ebessel function, eb_w from ebessel function, either 'expcos' or 'ebessel' depending on which is the best function);].
%           quadtree_n_points = String. Name of a text file with the number of pixels combined into each datapoint, if using quadtree. This is used to calculate the var-covar matrix. If not using quadtree this is a text file with '1's in each row, the same length as the InSAR datafile.
%           GPS_datafile = String. 'yes' or 'no'.
%           InSAR_coordinate_unit = 'utm' or 'long/lat'. Everything will be converted to UTM as necessary.
%           atolls_datafile = String. Name text file with uplift only data in.
%           atolls_coordinate_unit = String. 'utm' or 'long/lat'
%           UTMzone = Number (ignoring letter). The UTM zone e.g. 10
%           origin = A long/lat in the centre of your area [-124 38]. This is only used if your data spans multiple UTM zones and in that case everything is converted to a local coordinate system, using the origin.
%
%     'testing' is a structure, ONLY USED FOR TESTING. It allows you to add noise and load the true model to compare to your inversion solution:
%           testing_mode = String. 'yes' or 'no'. Note that if set to 'yes' any file in data.InSAR_datafile or data.GPS_datafile will be ignored and instead data in testing.making_model will be used.
%           making_model = String. File of saved data from 'making_model.m' 
%           var = Number. Imaginary values of variance for dataset. Assume no covariance. 
%           add_noise = Number. 0 = no noise, 1 = 1sigma noise, 2 = 2sigma-noise, etc etc
% 
%     'invert' is a structure with details of how you wish the inversion to be performed:
%           inversion_type = String. 'least_squares' (shame on you), or 'bayesian'
%           smoothing = String. 'None' gives no smoothing, 'laplacian' does Laplacian smoothing, 'VK' does von Karman smoothing.
%           smooth_across_fault_strands = String. 'yes' smoothes across different fault strands as if it were one fault. 'no' treats each fault as a different earthquakes. [] means I'll calculate which option you should choose based on the fault geometry .
%           iterations = Number. Number of iterations, only relevant for Bayesian inversions.
%           regularise_moment = String. 'yes' or 'no'. This adds an M0 prior likelihood - a normal distribution, with mean and std given in 'data' structure. 
%           solve_for_correlation_length = String. 'yes' or 'no'. This is the option to solve for the along-strike and down-dip correlation terms for VK smoothed solutions. This is not functional yet.
%           solve_for_dip = String. 'yes' or 'no'.
%           solve_for_fault_size = String. 'yes' or 'no'. NOTE: THIS ONLY WORKS WITH ONE FAULT STRAND, currently.
%           slip_initial = Number. Estimate of slip values - all slip patches will be assigned this amount of slip initially, with a bit of noise added.
%           step_size = Value. Maximum stepsize, plus or minus, allowed in your random walk when generating a slip_trial, for Bayesian inversions. This is used initially and then is adjusted through the iteration llh2local. For 'invert.slip_prior = boxcar', step_size is in metres. For 'invert.slip_perior = logarithmic' it's like a percentage of slip, or something.
%           alpha2_initial = Number. One number for each fault strand, of initial alpha^2 (variance) value. This is adjusted through the iteration.
%           alpha2_step_size = Number. Initial step size of alpha2_modelparameter, which is then adjusted through the inversion.
%           probability_target_initial = Number. Initial probability target - if we add half the step size to a parameter, this is the probability we're aiming the perturbation to make. This is adjusted through the iteration process: if the rejection rate is too high, then the probability decreases to try to decrease step sizes. % JUST FOR NOW while I'm sorting out sensitivity
%           pad_edge_with_zeros = String. 'yes' or 'no' - whether you'd like to impose zero slip around bottom, left and right of your predefined fault. This will affect the VK correlation most. NOT THOROUGHLY TESTED - probably won't work.
%           variable_rake = String. 'yes' or 'no'
%           load_old_MCMC_chain = String. Either a name of a saved file, or 'no'. If this is a name of a file, then the MCMC chain will continue from the past max-likelidhood solution, with same step sizes, alpha^2, and probability target. BE CAREFUL using this if you're loading a chain with DIFFERENT inversion parameters. some things may not work.
%           atolls_datafile: 'none'
%           solve_for_InSAR_offset: String. 'no' or 'yes'
%           solve_for_fault_size: String. 'yes' or 'no'
%           solve_for_beta = String. 'yes' or 'no' - beta is a hyperparameter on the data.
%           beta_initial = Number. Starting value of hyperparameter that acts on the var-covariance matrix. If this is set to 1 and beta_step_size is set to 0 then beta will stay at 1 throughout the whole inversion.
%           beta_step_size = Number. Starting value of hyperparmaeter stepsize. Default is 0.
%           simulated_annealing_start = String. 'yes' or 'no'. slipBERI will perform an initial simulated annealing inversion, to use as the starting parameter for the Bayesian inversion.
%           select_patches = String. This is the name of the matlab file that contains the numbers of which slip patches you wish to be 'on' during this inversion.
%           circular_harmonics = String. 'yes' or 'no'. If invert.solve_for_fault_size is 'yes' this solves for the size of the fault using circular harmonics.
%           ensemble_sampling = String. 'yes' or 'no' for whether you want to explode parameter space using ensemble sampling, rather than by adding a random number to the model parameters. NOT FULLY FUNCTIONAL.
%           ensemble_move_type = String. If using ensemble_sampling to explore parameter space then you must select 'stretch' or 'walk' for the method of doing so. NOT FULLY FUNCTIONAL.
%           quickcheck = String. 'yes' or 'no'. Before commencing the inversion, your data will be plotted for you to click on to say that it's okay.
%           add_correlation_matrix_stabiliser = String. 'yes' or 'no'. Add a small term (akin to a nugget) to the diagonal of the von Karman matrix sigma_s. Unlikely to be necessary unless solving for the fault size.
%
%     'priors' is a structure containing the priors
%           slip_prior = String. 'boxcar' or 'gaussian' or 'logarithmic', for Bayesian inversions. Note that if you use logarithmic, you will probably need to increase your number of iterations, since there are more rejections, and also decrease your step size (maybe 0.02 is reasonable?).
%           min_slip = Number. Minimum value of slip allowed, in metres. For Bayesian inversions. 
%           max_slip = Number. Maximum value of slip allowed, in metres. For Bayesian inversions. 
%           alpha2_prior = String. 'logarithmic' or [], depending on if you want to use a logarithmic prior on alpha2 or not.
%           alpha2_flag = String. Either 'bothsame' if you wish to use the same alpha2 hyperparameter on multiple fault strands or [] if not. 
%           max_alpha2 = Number. Maximum permitted alpha^2 hyperparameter. Note that if using logarithmic prior this is the maximum alpha^2 value, not maximum 10^(alpha^2) value.
%           min_alpha2 = Number. Minimum permitted alpha^2 hyperparameter. Note that if using logarithmic prior this is the minimum alpha^2 value, not maximum 10^(alpha^2) value.
%           min_rake = Number. Minimum value of permitted rake. Rake 0 = left-lateral strike-slip. Rake 180 = right-lateral strike-slip. Rake 90 = thrust. Rake -90 (or 270) = normal. Can be a matrix of one value per each fault_strand_for_smoothing, or if just one value used then same value is used across all fault strands.
%           max_rake = Number. Maximum value of permitted rake. Rake 0 = left-lateral strike-slip. Rake 180 = right-lateral strike-slip. Rake 90 = thrust. Rake -90 (or 270) = normal. Can be a matrix of one value per each fault_strand_for_smoothing, or if just one value used then same value is used across all fault strands.
%           min_dip = Number. Minimum value of dip permitted. If there are multiple fault strands, can enter one value per fault strand. If only one number is entered this is used for all fault strands.
%           max_dip = Number. Maximum value of dip permitted. If there are multiple fault strands, can enter one value per fault strand. If only one number is entered this is used for all fault strands.
%           predominant_faulting_style: String in curly brackets. Either {'ss'} or {'ds'}. This is used in von Karman smoothing to calculate the along-strike and down-dip correlation lengths.
%           max_offset = Number. Maximum permitted offset if invert.solve_forInSAR_offset = 'yes'.
%           min_offset = Number. Minimum permitted offset if invert.solve_forInSAR_offset = 'yes'.
%           min_beta = Number. Minimum value of permitted beta.           
%           max_beta = Number. Maximum value of permitted beta.
%           min_circharm_coeffs = Number. Minimum permitted circharm coefficient.
%           max_circharm_coeffs = Number. Maximum permitted circharm coefficient.
%           min_circharm_phi = Number. Minimum permitted circharm rotation (radians).
%           max_circharm_phi = Number. Maximum permitted circharm rotation (radians).
%           min_circharm_center = Number. Maximum permitted circharm center location (m).
%           max_circharm_center = Number. Maximum permitted circharm center location (m).
% 
%     'elastic_params' is a structure containing the Lame elastic constants whic are required for kernel calculation (Okada 1985)
%         lambda=3.23e10 ;
%         mu_okada=3.23e10 ;
% 
%     'display' is a structure which details how you would like your data to be displayed.
%         plotmean = String, 'yes' or 'no'. whether you'd like to plot the mean and standard deviation or not.
%         plotmode = String, 'yes' or 'no'. whether you'd like to plot the mode and standard deviation or not.
%         plotmedian = String, 'yes' or 'no'. whether you'd like to plot the median and standard deviation or not.
%         plotmaxlikely = String, 'yes' or 'no'. whether you'd like to plot the maximum likelihood and standard deviation or not.
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
%     'housekeeping' is important for keeping your files (as well as your mind, body and soul) tidy:
%           save_name = String. name you'd like to call your data run when you save. Note that slipBERI automatically appends it with your inversion style
%
% *********** % ************ % ********** % ********** % ********** % *****
%
% This function relies on many scripts, found in 'mbin'. Each script should
% give full details of what it requires, does, and outputs, found by typing
%                          'help <name_of_script>'
%
% Ruth Amey (rmja) 18th-Nov-2014 onwards...
% With much help from Andy Hooper & some scripts from GJ Funning, TJ Wright, D Bekaert and T Ingleby
%
% When using this code please cite:
% Amey, R. M. J., Hooper, A., & Walters, R. J. (2018). A Bayesian method 
% for incorporating self‐similarity into earthquake slip inversions. 
% Journal of Geophysical Research: Solid Earth, 123, 6052–6071. 
% or
% Amey, R.M.J., Hooper, A. and Morishita, Y. Going to Any Lengths: Solving
% for Fault Size and Fractal Slip for the 2016, Mw 6.2 Central Tottori
% Earthquake, Japan, using a Trans-dimensional Inversion Scheme, Journal of
% Geophysical Research: Solid Earth (IN REVIEW)


%% CHANGE HISTORY - for the backseat drivers amongst you ******************

% 'Date, initial, details' if you please. Let's keep it tidy in here.

% deets: rmja (eermja@leeds.ac.uk / R.M.J.Amey@leeds.ac.uk)

%  18-nov-2014  rmja  Started writing
%  26-feb-2015  rmja  Code actually ran the whole way through
%  08-may-2015  rmja  Moved from a script to a function
%  29-jul-2015  rmja  Fixed plotting InSAR, GPS and residuals
%  31-jul-2015  rmja  Added variable rake
%  18-aug-2015  rmja  Added converting GPS from lat/long to UTM
%  24-aug-2015  rmja  changed 'logprob_VK_temp' to 'logprior_VK_temp'; 'logprob_VK_trial' to 'logprior_VK_trial'; and 'logprob_VK_curr' to 'log_prior_VK_curr'.
%   7-sep-2015  rmja  Solving for alpha_squared as a hyperparameter
%  28-sep-2015  rmja  Fixed sensitivity - aiming for an ideal probability perturbation, which is also solved for
%   7-oct-2015  rmja  Add NaN checks
%   8-oct-2015  rmja  Added zero padding option on left, right, bottom edges of fault.
%  22-oct-2015  rmja  Corrected calculating RMS to WRMS
%  27-oct-2015  rmja  Separated 2D and 3D GPS
%   2-nov-2015  rmja  Can now solve for multiple fault strands
%   9-nov-2015  rmja  Remove points a mindist from the fault
%   5-jan-2016  rmja  Added Laplacian smoothing into MCMC. changed logprior_VK_curr to logprior_curr, logprior_VK_temp to logprior_temp and logprior_VK_trial to logprior_trial
%   5-jan-2016  rmja  Can choose whether to smooth across multiple fault strands or not
%  14-jan-2016  rmja  Added in M0 regularisation
%  10-feb-2016  rmja  Solve for along-strike and down-dip Vk correlation lengths
%  26-feb-2016  rmja  Solve for dip
%  16-jun-2016  rmja  Hacked process_faultdata_center a LOT so that when I plot results using imagesc it always plots the northerly and westerly patch in the top left.
%   4-jul-2016  rmja  Correcting weighting different datasets
%   8-jul-2016  rmja  Sorted out utm/latlong/local coordinates
%  13-jul-2016  rmja  Added in atolls as a dataset
%  14-jul-2016  rmja  Added in solving for InSAR offset
%   4-aug-2016  rmja  Fix solving for multiple faults under the surface (ramp flat ramp)
%   9-dec-2016  rmja  Added in solving for fault length
%   1-feb-2017  rmja  Changed varcovar calculation to calculating variogram
%  13-feb-2017  rmja  Fixed smoothing so you can identify which fault strands you want to smooth across and which you don't (whereas previously was smooth across all or smooth across none)
%  14-feb-2017  rmja  Fixed solving for patches being on/off so you can now use it with multiple fault strands
%     4-aug-2017  rmja  So you can choose different min/max rakes on different fault strands
%    23-aug-2017  rmja  Solve for fault size with circular harmonics
%    25-aug-2017  rmja  Ensemble approach to taking step sizes (Goodman and Weare, 2010) --- didn't really work for high number or model parameters
%    15-mar-2018  rmja  Parallelising ensemble move (Foreman-Mackey et al., 2013) --- didn't really work for high number or model parameters
%    27-mar-2018  rmja  Tidied code and fixed all the archaeic nonsense, lurking in subfunctions
%    22-may-2018  rmja  Solving for correlation lenghts a_as and a_dd from their distributions, when solving for fault size.

%% CHECK PATHS

% WRITE THIS
% *************************************************************************
profile on


%% Let's check our inputs are all fine-and-dandy before we proceed

%addpath(genpath('/nfs/see-fs-01_users/eermja/Documents/MATLAB/slipBERImbin')); 
run('error_messages.m');


%% FUN WITH FUNCTIONS - sit back, relax, let the magic begin **************

tic

% *********** % ************ % ********** % ********** % ********** % *****
% Introductions
% *********** % ************ % ********** % ********** % ********** % *****
disp(' ');
disp('Now running slipBERI');
disp('Written by Ruth Amey, with much help from Andy Hooper and others...');
disp(' ')

rng('shuffle')      % Fun fact - if you don't set your rng to shuffle then every time you open matlab fresh it'll give you the same random numbers in a row. Try it - I bet you £50 if you open a new matlab and type 'rand' you'll get 0.8147


% *********** % ************ % ********** % ********** % ********** % *****
% Look at the data
% *********** % ************ % ********** % ********** % ********** % *****
use_local_coordinate_system = 'no'; % by default use utmx utmy, not a local coordinate system

if strcmp(testing.testing_mode, 'yes') == 1
    load(testing.making_model, 'locs');
    n_datasets = 1;
    los_vector = ones(1,length(locs));
    data.GPS_datafile = 'none';
    first_InSAR_scene_numbers = [];
    InSAR_identifyer = [];
    n_InSAR_scenes = [];
else
    
    n_datasets = 0;
    
    % Sort out InSAR - data, locations, line of sight vector
    if strcmp(data.InSAR_datafile, 'none') ~= 1       % if it's TRUE that there is NOT 'none' InSAR data. i.e. if there IS an insar datafile
        n_datasets = n_datasets + 1;

        d_InSAR = [];
        locs_InSAR = [];
        los_vector_InSAR = [];
        counter = 1;
        InSAR_identifyer = [];
        n_InSAR_scenes = size(data.InSAR_datafile,2);
        first_InSAR_scene_numbers = 1;
        for p = 1:n_InSAR_scenes
            [ locs_tmp, d_tmp, los_vector_tmp ] = read_InSAR( cell2mat(data.InSAR_datafile(p)));
            d_InSAR = [d_InSAR, d_tmp];
            locs_InSAR = [locs_InSAR, locs_tmp];
            los_vector_InSAR = [los_vector_InSAR, los_vector_tmp];
            InSAR_identifyer = [InSAR_identifyer; counter*ones(length(d_tmp),1)];
            first_InSAR_scene_numbers = [first_InSAR_scene_numbers; (length(d_tmp)+first_InSAR_scene_numbers(end))];
            counter=counter+1;
            locs_tmp=[];
            los_vector_tmp=[];
            d_tmp=[];
        end
        first_InSAR_scene_numbers(end) = [];
        
        if strcmp(data.InSAR_coordinate_unit, 'long/lat') == 1
            locs_InSAR_latlong = locs_InSAR;
            [locs_InSAR(1,:),locs_InSAR(2,:)]=ll2utm(locs_InSAR(2,:),locs_InSAR(1,:));
        elseif  strcmp(data.InSAR_coordinate_unit,'utm') == 1
            locs_InSAR = locs_InSAR*1000;
            [locs_InSAR_latlong(2,:), locs_InSAR_latlong(1,:)] = utm2ll( locs_InSAR(1,:), locs_InSAR(2,:), repmat(data.UTMzone, 1, length(locs_InSAR)));
        end
        
        if locs_InSAR_latlong(1,:) < 0  % want long from
            locs_InSAR_latlong(1,:) = locs_InSAR_latlong(1,:) + 360;
        end

        clear d_tmp
        
    else
        d_InSAR = [];
        locs_InSAR = [];
        n_InSAR_scenes = [];
        first_InSAR_scene_numbers = [];
        InSAR_identifyer = [];
    end
     
    % Sort out GPS - data, locations, line of sight vector
    if strcmp(data.GPS_datafile_2d, 'none') + strcmp(data.GPS_datafile_3d, 'none') ~= 2        % if it's TRUE that there is NOT 'none' GPS data. i.e. if there IS an GPS datafile
        n_datasets = n_datasets + 1;        % GPS just counts as one dataset, even if 2d and 3d. yup.
        
        if strcmp(data.GPS_datafile_2d, 'none') ~= 1   % 2D GPS
            [ locs_GPS_2d, d_GPS_2d, los_vector_GPS_2d, sigma_gps_2d ] = read_GPS( data.GPS_datafile_2d );
            d_GPS_e_2d = d_GPS_2d(1:2:length(d_GPS_2d));
            d_GPS_n_2d = d_GPS_2d(2:2:length(d_GPS_2d));
            d_GPS_up_2d = [];
        else
            d_GPS_2d = [];
            d_GPS_e_2d = [];
            d_GPS_n_2d = [];
            d_GPS_up_2d = [];
            locs_GPS_2d = [];
            los_vector_GPS_2d = [];
            sigma_gps_2d = [];
        end
        if strcmp(data.GPS_datafile_3d, 'none') ~=1   % 3D GPS
            [ locs_GPS_3d, d_GPS_3d, los_vector_GPS_3d, sigma_gps_3d ] = read_GPS( data.GPS_datafile_3d );
            d_GPS_e_3d = d_GPS_3d(1:3:length(d_GPS_3d));
            d_GPS_n_3d = d_GPS_3d(2:3:length(d_GPS_3d));
            d_GPS_up_3d = d_GPS_3d(3:3:length(d_GPS_3d));
        else
            d_GPS_3d = [];
            d_GPS_e_3d = [];
            d_GPS_n_3d = [];
            d_GPS_up_3d = [];
            locs_GPS_3d = [];
            los_vector_GPS_3d = [];
            sigma_gps_3d = [];
        end

            d_GPS = [ d_GPS_2d, d_GPS_3d];
            d_GPS_e = [ d_GPS_e_2d, d_GPS_e_3d];
            d_GPS_n = [ d_GPS_n_2d, d_GPS_n_3d];
            d_GPS_up = [ d_GPS_up_2d, d_GPS_up_3d];
            locs_GPS = [ locs_GPS_2d,locs_GPS_3d ];
            los_vector_GPS = [los_vector_GPS_2d ,los_vector_GPS_3d ];
            sigma_gps = [ sigma_gps_2d ; sigma_gps_3d ];
            data.GPS_datafile = 'yes';
        
        if strcmp(data.GPS_coordinate_unit, 'long/lat') == 1
                [locs_GPS(1,:),locs_GPS(2,:), ~ ] = ll2utm(locs_GPS(2,:), locs_GPS(1,:));
                if  strcmp(data.GPS_datafile_2d, 'none') ~= 1
                    [locs_GPS_2d(1,:),locs_GPS_2d(2,:), ~ ] = ll2utm(locs_GPS_2d(2,:), locs_GPS_2d(1,:));
                end
                if  strcmp(data.GPS_datafile_3d, 'none') ~=1
                    [locs_GPS_3d(1,:),locs_GPS_3d(2,:), ~ ] = ll2utm(locs_GPS_3d(2,:), locs_GPS_3d(1,:));
                end
        else
            %locs_GPS = locs_GPS*1000;
        end
        
        locs_GPS_unique(1,:) = unique(locs_GPS(1,:), 'stable');
        locs_GPS_unique(2,:) = unique(locs_GPS(2,:), 'stable');
           
    else 
        data.GPS_datafile = 'none';
        d_GPS = [];
        locs_GPS = [];
        los_vector_GPS = [];
    end
    
    % Sort out atolls, if you have atolls
    if strcmp(data.atolls_datafile, 'none') ~= 1
        [ locs_atolls, d_atolls, los_vector_atolls, sigma_atolls ] = read_GPS( data.atolls_datafile );
        % Convert from latlong and name something different
        if strcmp(data.atolls_coordinate_unit,'long/lat') == 1
            locs_atolls_latlong = locs_atolls;
            [locs_atolls(1,:),locs_atolls(2,:), ~ ] = ll2utm(locs_atolls(2,:), locs_atolls(1,:));
        end
    else
        d_atolls = [];
        locs_atolls = [];
        los_vector_atolls = [];
    end
        


%         % see if you need to convert to a local coordinate system or keep as utm
%         
%          if strcmp(data.GPS_datafile, 'none') ~= 1 && strcmp(data.GPS_coordinate_unit, 'long/lat') ==1
%              [~,~, zones_GPS ] = ll2utm(locs_GPS(2,:), locs_GPS(1,:));
%          else
%              zones_GPS = [];
%          end
%         
%         if strcmp(data.InSAR_datafile, 'none') ~= 1
%             [~,~, zones_InSAR ] = ll2utm(locs_InSAR_latlong(2,:), locs_InSAR_latlong(1,:));
%         else
%             zones_InSAR = []; 
%         end
%         
%         if strcmp(data.atolls_datafile, 'none') ~= 1 
%             if strcmp(data.atolls_coordinate_unit, 'long/lat') == 1
%                 [~,~, zones_atolls ] = ll2utm(locs_atolls_latlong(2,:), locs_atolls_latlong(1,:));
%             else
%                 zones_atolls = []; 
%             end   
%         else
%             zones_atolls = [];
%         end
%         
%         zones_test = [zones_GPS, zones_InSAR, zones_atolls];
% 
%         if any(diff(zones_test) ~= 0 )
%             disp('having commented out converting to local coordinate system if UTM zones are different')
%             use_local_coordinate_system = 'yes';
%             disp('Converting GPS and InSAR to local coordinate sytem...')
%             
%             if strcmp(data.GPS_datafile, 'none') ~= 1 
%                 llh_GPS = [locs_GPS(1,:); locs_GPS(2,:); zeros(1,length(locs_GPS))];
%                 xy_GPS = llh2local(llh_GPS, data.origin);
%                 locs_GPS(1,:) = xy_GPS(1,:)*1000;   % convert to meters
%                 locs_GPS(2,:) = xy_GPS(2,:)*1000;   % convert to meters               
%             end
%             
%             if strcmp(data.InSAR_datafile, 'none') ~= 1
%                 llh_InSAR = [locs_InSAR(1,:); locs_InSAR(2,:); zeros(1,length(locs_InSAR(1,:)))];
%                 xy_InSAR = llh2local(llh_InSAR, data.origin);
%                 locs_InSAR(1,:) = xy_InSAR(1,:)*1000;   % convert to meters
%                 locs_InSAR(2,:) = xy_InSAR(2,:)*1000;   % convert to meters
%             end
%       
%         end
        
    % Combine all data, locations, line of sights
    d = [d_InSAR, d_GPS, d_atolls];
    locs = [locs_InSAR, locs_GPS, locs_atolls];
    los_vector =  [los_vector_InSAR, los_vector_GPS, los_vector_atolls];
    
    d = d';
    
    n_data = length(d);
    d_identifyer = [ 1* ones(length(d_InSAR),1); 2*ones(length(d_GPS),1); 3*ones(length(d_atolls),1)];


end


% *********** % ************ % ********** % ********** % ********** % *****
% Sort out the fault plane
% *********** % ************ % ********** % ********** % ********** % *****

% Create dislocation model for each patch
[disloc_model, spatial_model1, spatial_model2, spatial_model3, spatial_model4, fault_coords, ~, n_along_strike_patches, n_down_dip_patches, fault_length, fault_width,n_fault_strands, strike, fault_strand_togetherness ] =...
	process_faultdata_centre(fault, invert, use_local_coordinate_system, data.origin, testing);             % if you're padding the edges with zeros, THESE EXTRA PATCHES GET ADDED HERE

spatial_model2column = reshape(spatial_model2, [],1);
spatial_model3column = reshape(spatial_model3, [],1);

% %Save for GMT
% fault_coords_latlong = [];
% for i = 1:n_fault_strands
%     [LAT,LON]=utm2ll([fault_coords(i,1), fault_coords(i,3)],[fault_coords(i,2), fault_coords(i,4)],data.UTMzone); 
%     fault_coords_latlong = [fault_coords_latlong; (LON)', LAT'];
% end
% dlmwrite('fault_coords_latlong.gmt', fault_coords_latlong);

% Start with some easy thing
total_n_slip_patches = size(disloc_model, 2);
n_slip_patches_on_each_fault_strand = n_along_strike_patches .* n_down_dip_patches;      % one value per fault strand
first_patch_in_strand(1) = 0;
first_patch_in_strand(2:(n_fault_strands+1),1) = cumsum(n_slip_patches_on_each_fault_strand);
first_patch_in_strand = first_patch_in_strand+1;
first_patch_in_strand(n_fault_strands+1) = [];                          % so that you know which number slip patch is the first one in each new strand. comes in handy
last_patch_in_strand = cumsum(n_slip_patches_on_each_fault_strand);     % so that you know which number slip patch is the last one in each new strand. comes in handy
onoffidentifyer = ones(total_n_slip_patches,1);

along_strike_sep_dist = zeros(n_fault_strands,1);
down_dip_sep_dist = zeros(n_fault_strands,1);
fault_strand_identifyer = [];
fault_strand_identifyer_for_smoothing = [];
difference_between_faultssegments_and_smoothingfaultsegments = diff(fault_strand_togetherness);
difference_between_faultssegments_and_smoothingfaultsegments = [ 0;difference_between_faultssegments_and_smoothingfaultsegments];
counter = 1;
for j = 1 : n_fault_strands
    along_strike_sep_dist(j,1) = spatial_model2(first_patch_in_strand(j)); % one number for each strand. ASSUMING ALL FAULT PATCHES ARE THE SAME SIZE IN ONE FAULT STRAND. otherwise I can make a matrix of separation distances and separate terms for each fault patch
    down_dip_sep_dist(j,1) = spatial_model3(first_patch_in_strand(j));     % one number for each strand. ASSUMING ALL FAULT PATCHES ARE THE SAME SIZE IN ONE FAULT STRAND.
    identifyer = j * ones( n_slip_patches_on_each_fault_strand(j),1);
    secondidentifyer = counter * ones( n_slip_patches_on_each_fault_strand(j),1);
    fault_strand_identifyer = [fault_strand_identifyer; identifyer ];
    if difference_between_faultssegments_and_smoothingfaultsegments(j)==1      % then this is considered a new strand to the last strand
        counter = counter + 1;
        secondidentifyer = counter * ones( n_slip_patches_on_each_fault_strand(j),1);
        fault_strand_identifyer_for_smoothing = [fault_strand_identifyer_for_smoothing; secondidentifyer ];     % one value for each fault patch
    elseif difference_between_faultssegments_and_smoothingfaultsegments(j)==0  % then this is considered the same strand as the last strand
        fault_strand_identifyer_for_smoothing = [fault_strand_identifyer_for_smoothing; secondidentifyer ];     % one value for each fault patch
    end
 
end

clear identifyer
    
% If padding the edges with zero, put zero slip on these patches.
if strcmp(invert.pad_edges_with_zeros, 'yes') ==1
    disloc_model(6, 1: n_down_dip_patches) = 0;                             % add zeros on left of fault.
    disloc_model(6, (total_n_slip_patches-n_down_dip_patches):end) = 0;     % add zeros on right of fault.
    disloc_model(6, n_down_dip_patches:n_down_dip_patches:end) = 0;         % add zeros on bottom of fault.
end

if any(disloc_model(8,first_patch_in_strand)~=0)    % if any of the faults start below ground then we need to sum all the down dip patches, and they'll all have the same number of along strike patches
    total_n_along_strike_patches = n_along_strike_patches(1);       % this is necessary because sometimes we have multiple strands laterally (e.g. napa) and sometimes we have them going down depth (e.g. nepal)
    total_n_down_dip_patches = sum(n_down_dip_patches);
else    % if all the faults start at the surface but there are a few extending laterally then we need to sum all the alongstrike patches, and they'll all have the same number of downdip patches
    total_n_along_strike_patches = sum(n_along_strike_patches);
    total_n_down_dip_patches = n_down_dip_patches(1);
end

     
% Multiple fault strands...............................................

    if strcmp(invert.smooth_across_fault_strands, []) ==1 % if they left it to me to decide whether to treat fault strands differently or not.       
        if any(abs(diff(strike))) > 18                             % austin elliott et al 2015 - bends greater than 18 degrees terminate rupture, so should treat as different faults
            invert.smooth_across_fault_strands = 'yes';
        else
            invert.smooth_across_fault_strands = 'no';
        end
    end
    
    if strcmp(invert.smooth_across_fault_strands, 'yes') == 1       % note that if you do this then only the fault strands with the SAME IDENTIFYER IN the fault geometry file will be smoothed across
        %n_fault_strands_for_smoothing = 1;
        n_fault_strands_for_smoothing = length(unique(fault_strand_togetherness));
        %fault_strand_identifyer_for_smoothing = []; % we made this assuming all the different fault strands would be separate. so we'll ignore that now
        
        difference_between_faultssegments_and_smoothingfaultsegments = diff(fault_strand_togetherness);
        
        for r = 1: n_fault_strands_for_smoothing  % this works out how many fault strands you wanna have
            fault_length_for_smoothing(r,1) = sum(fault_length(fault_strand_togetherness==r)); % I want this to sum faults that are going to be considered together
            %fault_width_for_smoothing(r,1) = sum(fault_width(fault_strand_togetherness==r));
            fault_width_for_smoothing(r,1) = fault_width(r,1); % assuming that all faults are the same depth
            first_patch_in_strand_for_smoothing_master(r,1) = find(fault_strand_identifyer_for_smoothing==r,1);
            last_patch_in_strand_for_smoothing_master(r,1) = find(fault_strand_identifyer_for_smoothing==r,1, 'last');
            n_down_dip_patches_for_smoothing(r,1) = total_n_down_dip_patches; % because we're treating it as one fault. assuming all faults have some number of down dip patches
            n_along_strike_patches_for_smoothing(r,1) = sum(n_along_strike_patches(fault_strand_togetherness==r)); % because we're treating it as one fault
            n_slip_patches_on_each_fault_strand_for_smoothing(r,1) = sum(fault_strand_identifyer_for_smoothing==r);      
        end
    elseif strcmp(invert.smooth_across_fault_strands, 'no') == 1
        n_fault_strands_for_smoothing = n_fault_strands; 
        n_down_dip_patches_for_smoothing = n_down_dip_patches; % because we're treating it as separate faults as described in fault descriptor file
        n_along_strike_patches_for_smoothing = n_along_strike_patches; % because we're treating it as separate faults as described in fault descriptor file  
        first_patch_in_strand_for_smoothing_master = first_patch_in_strand;
        last_patch_in_strand_for_smoothing_master = last_patch_in_strand;
        fault_length_for_smoothing = fault_length; % this is one value for each separate strand, regardless of whether attached or not
        fault_width_for_smoothing = fault_width; % this is one value for each separate strand, regardless of whether attached or not
        n_slip_patches_on_each_fault_strand_for_smoothing = n_slip_patches_on_each_fault_strand;
    end
 
first_patch_in_strand_for_smoothing = first_patch_in_strand_for_smoothing_master;
last_patch_in_strand_for_smoothing = last_patch_in_strand_for_smoothing_master;
n_slip_patches_ON_on_each_fault_strand_for_smoothing = n_slip_patches_on_each_fault_strand_for_smoothing; % because they're all on to start off with
previous_n_slip_patches_ON_on_each_fault_strand_for_smoothing = n_slip_patches_ON_on_each_fault_strand_for_smoothing;


   
% *********** % ************ % ********** % ********** % ********** % *****
% If in testing mode...
% *********** % ************ % ********** % ********** % ********** % *****
% Sort out synthetic data................................... 
if strcmp(testing.testing_mode, 'yes') == 1   
    % synthetic data.......................................................
        load(testing.making_model, 'u');        % these are the true observations, from the forward model of the true slip 
            d_E = u(1,:); 
            d_N = u(2,:);
            d_up = u(3,:);
            d = [d_E'; d_N'; d_up'];            % THIS IS d in E, N, up , as opposed to LOS
            n_data = length(d);
            d_identifyer = zeros(length(d), 1);
       
    % add noise, if you like...............................................
        if testing.add_noise ~= 0              % add_noise tells you the std of noise to add. add_noise == 1 means add 1 sigma noise, add_noise == 2 means add 2 sigma noise
            % calculate standard deviation
            sigmas = std(u');
            sigma_E = sigmas(1);
            sigma_N = sigmas(2);
            sigma_up = sigmas(3);
            % create an array of random numbers, in the sigma range. multiply by add_noise, so that you add the right number of sigmas
            noise_E = (-sigma_E + (sigma_E - -sigma_E) * rand( length(u),1)) * testing.add_noise;% create an array of random numbers, in the sigma range
            noise_N = (-sigma_N + (sigma_N - -sigma_N) * rand( length(u),1)) * testing.add_noise;% create an array of random numbers, in the sigma range
            noise_up = (-sigma_up + (sigma_up - -sigma_up) * rand( length(u),1)) * testing.add_noise;% create an array of random numbers, in the sigma range
            % update data vector
            d_E = u(1,:) + noise_E';
            d_N = u(2,:) + noise_N';
            d_up = u(3,:) + noise_up';
            d = [d_E'; d_N'; d_up']; 
        end   
        
        d_E = [];       % get rid of it so you don't get confused / mess it up with removing point a mindist from the fault
        d_N = [];      
        d_up = []; 
        
    % find out true moment.................................................
    
    true1 = load(testing.making_model, 'spatial_model2');
    true2 = load(testing.making_model, 'spatial_model3');
    true3 = load(testing.making_model, 'synthetic_slip');
    true3.synthetic_slip = reshape(true3.synthetic_slip, length(true3.synthetic_slip), 1);
    true4 = load(testing.making_model, 'total_n_slip_patches');
    
    M0_true = sum((elastic_params.mu_okada) * (reshape(true1.spatial_model2, true4.total_n_slip_patches,1) .* reshape(true2.spatial_model3, true4.total_n_slip_patches,1)) .* (true3.synthetic_slip));
    data.seismic_moment = M0_true;     
end
      



% *********** % ************ % ********** % ********** % ********** % *****
% Create the G kernel............................................
% *********** % ************ % ********** % ********** % ********** % *****
disp('... and calculating G, using Okada 1985...')

if strcmp(testing.testing_mode, 'yes') == 1    % OR IF ONLY USING GPS DATA, basically this is just not projected into LOS
    for i = 1:total_n_slip_patches
        [u_for_unit_slip, ~]=disloc3d3(disloc_model(:,i), locs, elastic_params.lambda,elastic_params.mu_okada) ;
        G_E(:,i) = u_for_unit_slip(1,:);
        G_N(:,i) = u_for_unit_slip(2,:);
        G_up(:,i) = u_for_unit_slip(3,:);
        G = [G_E; G_N; G_up];
    end
end

if strcmp(invert.solve_for_dip, 'yes') == 1
    dip_LUT = priors.min_dip : priors.max_dip;      % Make a dip look up table
elseif  strcmp(invert.solve_for_dip, 'no') == 1  
    dip_LUT = disloc_model(4);   
end
    
  
if strcmp(invert.variable_rake, 'yes') == 1
        G_ss = zeros( n_data, total_n_slip_patches, length(dip_LUT));
        G_ds = zeros( n_data, total_n_slip_patches, length(dip_LUT));
        for i = 1: length(dip_LUT)
            % Calculate G for purely strike-slip
            disloc_model_ss = disloc_model;
            disloc_model_ss(5,:) = 0;                 % set rake of every patch to 0, which is purely LEFT LATERAL STRIKE SLIP
            disloc_model_ss(4,:) = dip_LUT(i);          % look up for each value of dip in LUT
            %[ G_ss(:,:,i), ~, ~, ~ ] = calculate_G( total_n_slip_patches, testing, data, disloc_model_ss, locs_struct, elastic_params, los_vector_struct );
            for j = 1: total_n_slip_patches
                [u_for_unit_slip, ~]=disloc3d3(disloc_model_ss(:,j), locs, elastic_params.lambda, elastic_params.mu_okada) ;
                if strcmp(testing.testing_mode, 'yes') ==1
                    G_E_ss(:,j) = u_for_unit_slip(1,:);
                    G_N_ss(:,j) = u_for_unit_slip(2,:);
                    G_up_ss(:,j) = u_for_unit_slip(3,:);
                    G_ss = [G_E_ss; G_N_ss; G_up_ss];
                else
                    G_ss(:,j) = sum(u_for_unit_slip.*los_vector);
                end
            end
            
            % Calculate G for purely dip-slip
            disloc_model_ds = disloc_model;
            disloc_model_ds(5,:) = 90;                 % set rake of every patch to 90, which is purely THRUST
            disloc_model_ds(4,:) = dip_LUT(i);
            %[ G_ds(:,:,i), ~, ~, ~] = calculate_G( total_n_slip_patches, testing, data, disloc_model_ds, locs_struct, elastic_params, los_vector_struct );
            for j = 1: total_n_slip_patches
                [u_for_unit_slip, ~]=disloc3d3(disloc_model_ds(:,j), locs, elastic_params.lambda, elastic_params.mu_okada) ;
                if strcmp(testing.testing_mode, 'yes') ==1
                    G_E_ds(:,j) = u_for_unit_slip(1,:);
                    G_N_ds(:,j) = u_for_unit_slip(2,:);
                    G_up_ds(:,j) = u_for_unit_slip(3,:);
                    G_ds = [G_E_ds; G_N_ds; G_up_ds];
                else
                    G_ds(:,j) = sum(u_for_unit_slip.*los_vector);
                end
            end
        end
else
    %[ G, G_E, G_N, G_up ] = calculate_G( total_n_slip_patches, testing, data, disloc_model, locs_struct, elastic_params, los_vector_struct );
    for i = 1: total_n_slip_patches
        [u_for_unit_slip, ~]=disloc3d3(disloc_model(:,i), locs, elastic_params.lambda, elastic_params.mu_okada) ;
        G(:,i) = sum(u_for_unit_slip.*los_vector);
    end
    G_ss = [];
    G_ds = [];    
end

 




% *********** % ************ % ********** % ********** % ********** % *****
% Variance-covariance matrix 
% *********** % ************ % ********** % ********** % ********** % *****
    % Set up data hyperparameter
    beta_step_size_initial = [];
    beta_ratio = 1; 
    
    if strcmp(invert.solve_for_beta, 'yes') == 1
        beta_initial = invert.beta_initial;
        beta_step_size_initial = invert.beta_step_size;
        beta_keep = zeros(1,ceil(invert.iterations));
        max_beta = priors.max_beta;
        min_beta = priors.min_beta;
    else
        max_beta = [];
        min_beta = [];
    end

if strcmp(testing.testing_mode, 'no') == 1
    
    %if you have InSAR data ..................................................
    if strcmp(data.InSAR_datafile, 'none') ~= 1
         disp('  ');
         disp(['Calculating var-covar matrix for Insar file ', data.InSAR_datafile]);

            % Assuming previously calculated sill and nugget and range, so load it.
            sigma_d_InSAR = [];
            for p = 1:size(data.varcovar_deets,2)
                [sill, nugget, variogram_range]=textread(cell2mat(data.varcovar_deets(p)),'%f %f %f');  % range in meters
                [npoints]=textread(cell2mat(data.quadtree_n_points(p)),'%f');  % range in meters
            
                % Calculate distance between all points
                x = locs_InSAR(1,InSAR_identifyer==p);   % meters
                y = locs_InSAR(2,InSAR_identifyer==p);
                x_dists =  pdist2(x',x');
                y_dists =  pdist2(y',y');
                dist = sqrt( x_dists.^2 + y_dists.^2);     % this is in meters

                % Calculate var-covar matrix
                sigma_d_tmp = (sill-nugget) * exp((-3*dist)/variogram_range)+diag(nugget./ npoints);      %if using quadtree.
                sigma_d_InSAR = blkdiag(sigma_d_InSAR, sigma_d_tmp);
            end
            
            % If you you want to weight it ....................................................
            sigma_d_InSAR = sigma_d_InSAR/data.weight_InSAR;
    else
         sigma_d_InSAR = [];
    end

    % if you have GPS data and you want to weight it  ....................................................
    if strcmp(data.GPS_datafile, 'none') ~= 1
      sigma_d_GPS = diag(sigma_gps.^2) / data.weight_GPS;   % variance is standard deviation squared. so if InSAR weight = 1 and GPS weight = 2, then we half the errors on GPS, coz we trust it quite as much
    else
        sigma_d_GPS = [];
    end
    
     % if you have atolls data and you want to weight it  ................................................
    if strcmp(data.atolls_datafile, 'none') ~= 1
      sigma_d_atolls = diag(sigma_atolls.^2) / data.weight_atolls;   % so if InSAR weight = 1 and GPS weight = 2, then we half the errors on GPS, coz we trust it quite as much
    else
        sigma_d_atolls = [];
    end
    
    sigma_d = blkdiag(sigma_d_InSAR, sigma_d_GPS, sigma_d_atolls);

elseif strcmp(testing.testing_mode, 'yes') == 1
    datastd = testing.std * ones(length(d),1);
    sigma_d = diag(datastd).^2;            % assume no covar, only diagonals
        
end

inv_sigma_d = inv(sigma_d);




 
%% Invert for slip - it's why we're all here, this is the good bit ********        

disp('     ');
disp('Solving for slip...');
disp('Oh boy, this is going to be fun.');
disp('    ');


%% OPTION 1 - least-squares inversion *************************************

% does least squares find the actual solution? a lesson on non-uniqueness
% no noise - least_squares finds correct solution. good.
% with noise - least-squares solution does not
% 
% if strcmp(invert.inversion_type, 'least-squares') == 1;
%    error('least-squares inversion isn''t actually an option'); 
% end

% AT THE MOMENT this doesn't solve for variable rake... G is calculated
% just using the rake given in fault.fault_descriptor_file

if strcmp(invert.inversion_type, 'least_squares') ==1 && strcmp(invert.smoothing, 'laplacian') == 1   
        disp('Finding Laplacian matrix, with smoothing factor = ') 
        [smooth_los_kernel,smooth_obsdispl_los] = add_smoothing(G, d, spatial_model1, spatial_model2, spatial_model3, spatial_model4, (invert.scalar_smoothing_factor^2), n_datasets) ; % borrowed from gjf. This function takes an already existing kernel and observation set, and adds a Laplacian-minimising smoothing operation to it.
    	G=smooth_los_kernel ;
   	 	d=smooth_obsdispl_los ;
end
            

if strcmp(invert.inversion_type, 'least_squares') == 1 && strcmp(invert.smoothing, 'tikhonov') +  strcmp(invert.smoothing, 'cv') ~= 1 
        disp('Solving for slip using fast non-negative least-squares...');
		ATA= G'* G;
    	ATD= G'*d ;
    	[least_sq_solution, ~]=fnnls(ATA,ATD) ;
        least_sq_residuals =  calc_WRMS_offset( least_sq_solution, G, d, inv_sigma_d, offset_curr );  
end
     
   
            

% For Tikhonov solution we want to minimise
%               (d-Gm)^2 + mu*(m)^2
% Which is the damped least squares solution with the special case of
% minimising the weighted misfit of hte data misfit and the model norm. This
% is the Tikhonov solution. To find the value of mu we solve the solution
% repeatedly till we find the value such that
%               (d-Gm)^2 = N * sigma^2
%               (d-Gm)   = sqrt(N) * sigma
% Where N is the number of degrees of freedom and sigma^2 is the variance
% (Since on average each error should equal sigma^2. So for N measurements
% the error should be N*sigma^2.)
% The solution for the damped least-squares solution is:
%           m_tik = inv(G'*G + mu*I) * (G'*d)
% Or if you want to include data error:
%    m_tik = inv(G' * inv_sigma_d * G + mu*I) * (G' * inv_sigma_d * d)


% if strcmp(invert.inversion_type, 'least_squares') == 1 && strcmp(invert.smoothing, 'tikhonov') == 1 
%         disp('Solving for slip using a linear least-squares tikhonov inversion...');
%         %G_tik = [G, ones(length(d),1) ]; % coz I want to solve for offset too
%         G_tik = G; % coz I want to solve for offset too
%         mew_tests = [0.1:0.1:10, 10:1:100];
%         tikhonov_criteria = total_n_slip_patches   * sigma_d(1); 
%         %tikhonov_criteria = sqrt(total_n_slip_patches)   * sqrt(sigma_d(1));
%         
%         for i= 1:length(mew_tests)
%             mew = mew_tests(i);
%             m_tik(:,i) = (G_tik' * inv_sigma_d * G_tik + mew * eye(total_n_slip_patches,total_n_slip_patches)) \ G_tik' * inv_sigma_d * d;      % since inv(G' * G + mu * I)*(G' * d) is slower and less accurate than (G' * G + mu * I)\(G' * d) 
%             resid(:,i) = sum((d - (G_tik*m_tik(:,i))).^2);
%             %resid(:,i) = sum(d - (G_tik*m_tik(:,i)));
%             tik_test(1,i) = resid(:,i)-tikhonov_criteria;
%         end
%         
%     % The best result is the one where (d-Gm) is closest to sqrt(N)*sigma
%     [~,I] = min(abs(resid-tikhonov_criteria));
%     best_mew = mew_tests(I);
%     best_m_tik = m_tik(:,I);
%     best_slip_tik = best_m_tik(1:end);
% 
%     figure; plot( mew_tests, resid )
%     ylabel('sum( d - Gm)')
%     xlabel('mew')
%     
%     best_offset_tik = 0;
% end 


%%%%%%%%%%%%%%% This is solving for rake as well %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% but no positivity constraint %%%%%%%%%%%%%%%%%%%%%%%%

% if strcmp(invert.inversion_type, 'least_squares') == 1 && strcmp(invert.smoothing, 'tikhonov') == 1
%         disp('Solving for slip and rake using a linear least-squares tikhonov inversion...');
%         G_tik = [G_ss, G_ds ]; 
%         mew_tests = [0:10:1000];
%         m_tik = zeros(total_n_slip_patches*2, length(mew_tests));
%         resid = zeros(1, length(mew_tests));
%         for i= 1:length(mew_tests)
%             mew = mew_tests(i);
%             m_tik(:,i) = (G_tik' * inv_sigma_d * G_tik + mew * eye(total_n_slip_patches*2,total_n_slip_patches*2)) \ (G_tik' * inv_sigma_d * d);      % since inv(G' * G + mu * I)*(G' * d) is slower and less accurate than (G' * G + mu * I)\(G' * d)
%             slip_tik = m_tik(:,i);
%             %offset_tik = m_tik(end);
%             least_sq_residuals =  calc_WRMS_offset( slip_tik, G_tik, d, inv_sigma_d, zeros(length(d),1) );  
%             resid(:,i) = sum((d - G_tik*m_tik(:,i)).^2);
%             %figure; faults= disloc_model; faults(6,:) = sqrt(m_tik(1:total_n_slip_patches,i).^2 + m_tik(total_n_slip_patches+1:end,i).^2); doplot3d(faults', 'jet');
%         end
%         
%     %The best result is the one where (d-Gm) is closest to sqrt(N)*sigma
%     tikhonov_criteria = (length(m_tik(:,1)))   * sigma_d(1);  
%     [~,I] = min(abs(resid-tikhonov_criteria));
%     best_mew = mew_tests(I);
%     best_m_tik = m_tik(:,I);
%     best_ss_slip_tik = best_m_tik(1:total_n_slip_patches);
%     best_ds_slip_tik = best_m_tik(total_n_slip_patches+1:end);
%     best_slip_tik = sqrt(best_ss_slip_tik.^2 + best_ds_slip_tik.^2);
%     best_rake_tik = 180 + atand( best_ds_slip_tik./best_ss_slip_tik);
%     
%            theta = 180 - best_rake_tik;
%            theta = degtorad(theta);
%            for r = 1 : total_n_slip_patches                                 % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake.
%                 best_G_tik(:,r) = cos(theta(r)) * G_ss(:,r) + sin(theta(r)) * G_ds(:,r);    % I THINK. check yo' trigonometry  
%            end
% 
%            best_offset_tik = 0;
% end

%%%%%%%%%%%%%%% This is tikhonov using fnnls %%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(invert.inversion_type, 'least_squares') == 1 && strcmp(invert.smoothing, 'tikhonov') == 1
        disp('Solving for slip and rake using a linear least-squares tikhonov inversion...');

        % We calculate Gss and Gds above, but I've copied this here, to
        % solve for 45 degrees around a specified rake
        centrerake = disloc_model(5,1);
        minrake = centrerake-45;
        maxrake = centrerake+45;
        
           disloc_model_minrake = disloc_model;
           disloc_model_minrake(5,:) = minrake;                 % set rake of every patch to 0, which is purely RIGHT LATERAL STRIKE SLIP
           disloc_model_minrake(4,:) = dip_LUT(i);          % look up for each value of dip in LUT
           [ G_minrake(:,:,i), ~, ~, ~ ] = calculate_G( total_n_slip_patches, testing, data, disloc_model_minrake, locs_struct, elastic_params, los_vector_struct );  
           % [ G_ss(first_patch_in_strand(i):last_patch_in_strand(i),first_patch_in_strand(i):last_patch_in_strand(i)), ~, ~, ~ ] = calculate_G( total_n_slip_patches, testing, data, disloc_model_ss, locs_struct, elastic_params, los_vector_struct );

           % Calculate G for purely dip-slip 
           disloc_model_maxrake = disloc_model;
           disloc_model_maxrake(5,:) = maxrake;                 % set rake of every patch to 90, which is purely THRUST
           disloc_model_maxrake(4,:) = dip_LUT(i);
           [ G_maxrake(:,:,i), ~, ~, ~] = calculate_G( total_n_slip_patches, testing, data, disloc_model_maxrake, locs_struct, elastic_params, los_vector_struct );

        
        mew_tests = 1:10:100;
        m_tik = zeros(total_n_slip_patches*2, length(mew_tests));
        resid = zeros(1, length(mew_tests));
        d_tik = [d; zeros((total_n_slip_patches*2),1)];
           
        weights = chol(inv(sigma_d));

        for i= 1:length(mew_tests)
            % Make new G using smoothing factor
            mew = mew_tests(i);
            G = [G_minrake, G_maxrake ; mew * eye(total_n_slip_patches*2) ];
            
            % Add weights and apply
            allW_tik = blkdiag(weights,eye(total_n_slip_patches*2)); % add in ones along diagonal for minimum norm
            G_tik = allW_tik*G;
            obs_tik = allW_tik*d_tik;

            % Solve for slip
            ATA=G_tik'*G_tik ;
            ATD=G_tik'*obs_tik ;
    		[m_tik(:,i), ~]=fnnls(ATA,ATD) ;

            least_sq_residuals =  calc_WRMS_offset( m_tik(:,i), [G_minrake, G_maxrake], d, inv_sigma_d, zeros(length(d),1) );  
            %resid(:,i) = sum((d - G_tik*m_tik(:,i)).^2);
            %figure; faults= disloc_model; faults(6,:) = sqrt(m_tik(1:total_n_slip_patches,i).^2 + m_tik(total_n_slip_patches+1:end,i).^2); doplot3d(faults', 'jet');
        end
        
    %The best result is the one where (d-Gm) is closest to sqrt(N)*sigma
    tikhonov_criteria = (length(m_tik(:,1)))   * sigma_d(1);  
    [~,I] = min(abs(least_sq_residuals-tikhonov_criteria));
    best_mew = mew_tests(I);
    best_m_tik = m_tik(:,I);
    best_minrake_slip_tik = best_m_tik(1:total_n_slip_patches);
    best_maxrake_slip_tik = best_m_tik(total_n_slip_patches+1:end);
    best_slip_tik = sqrt(best_minrake_slip_tik.^2 + best_maxrake_slip_tik.^2);
    best_rake_tik = minrake + atand( best_maxrake_slip_tik./best_minrake_slip_tik);
    best_rake_tik(best_minrake_slip_tik==0 & best_maxrake_slip_tik==0) = centrerake;
    
%            theta = 180 - best_rake_tik;
%            theta = degtorad(theta);
%            for r = 1 : total_n_slip_patches                                 % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake.
%                 best_G_tik(:,r) = cos(theta(r)) * G_ss(:,r) + sin(theta(r)) * G_ds(:,r);    % I THINK. check yo' trigonometry  
%            end
           best_G_tik = G_ss*diag(cosd(best_rake_tik)) +  G_ds*diag(sind(best_rake_tik));  
           

           best_offset_tik = 0;
end 





%%%%%%%%%%%%%%%% SVD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explanation of SVD goes here

if strcmp(invert.inversion_type, 'SVD') == 1 == 1
        disp('Solving for slip and rake using SVD...');
        G_svd = [G_ss, G_ds ];
        [U, S, V] = svd(G_svd);
        p = sum(diag(S)<1e-6);
        Sp = S(1:p, 1:p);
        Vp = V(:, 1:p);
        Up = U(:, 1:p);
        Gg = Vp * inv(Sp) * Up.';
            m_svd = Gg * d; 
            least_sq_residuals =  calc_WRMS_offset( m_svd, G_svd, d, inv_sigma_d, zeros(length(d),1) );  
        
           ss_slip_svd = m_svd(1:total_n_slip_patches);
           ds_slip_svd = m_svd(total_n_slip_patches+1:end);
           slip_svd = sqrt(ss_slip_svd.^2 + ds_slip_svd.^2);
           rake_svd = 180 + atand(ds_slip_svd./ss_slip_svd);
           
%            theta = 180 - rake_svd;
%            theta = degtorad(theta);
%            for r = 1 : total_n_slip_patches                                 % we wanna select the correct COLUMN, relating to the current slip patch, out of the G matrix for that value of rake.
%                 best_G_svd(:,r) = cos(theta(r)) * G_ss(:,r) + sin(theta(r)) * G_ds(:,r);    % I THINK. check yo' trigonometry  
%            end
           best_G_svd = G_ss*diag(cosd(rake_svd)) +  G_ds*diag(sind(rake_svd));  
           

           best_offset_tik = 0;
end 



%%%%%%%%%%%%%%%%% Cross validation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For cross-validation we pick which is the best smoothing factor by
% measuring the fit to missing data. A data point is dropped, the minimum
% norm solution is calculated for a given range of smoothing factors, and
% then the best smoothing factor is the one that predicts the missing data
% best.
% Tom has altered this (and sped it up) by instead of dropping one data
% point each time, he instead randomly reorders the data vector
% and cuts into k subsets. He uses k-1 subsets to estimate the inversion
% parameters and predict the data values at the remaining subset. This is
% done k times for each smoothing factor and the total misfit calculated by
% summing the residuals in each subset. The smoothing factor with the
% lowest overall misfit is the best option.


if strcmp(invert.inversion_type, 'least_squares') == 1 && strcmp(invert.smoothing, 'cv') == 1   
        
        % Use Tom's cross-validation code to select best smoothing factor
        G_cv = [G_ss, G_ds];
        smoothing_factors = [10^-3,10^-2,10^-1,10^0,10^1,10^2,10^3];
        k_subsets = 10;
    
        %[V_smoothing_factor,best_sf,figure_counter] = slipmodel_cross_validation(G,d,vcm,L,n_patches,smoothing_factors,k_subsets,figure_counter,project_str)
        [~,best_sf,~] = slipmodel_cross_validation_ruthhack(G_cv, d, sigma_d, total_n_slip_patches*2, smoothing_factors, k_subsets, 1, []);    % hacked from tom's laplacian version  

        disp(['Best smoothing factor is ', num2str(best_sf)])
        
        best_slip_cv = (G_cv' * sigma_d * G_cv + best_sf * eye(total_n_slip_patches*2,total_n_slip_patches*2)) \ G_cv' * sigma_d * d;
        best_slip_ss = best_slip_cv(1:total_n_slip_patches);
        best_slip_ds = best_slip_cv(total_n_slip_patches+1:end);
        best_slip_cv = sqrt(best_slip_ss.^2 + best_slip_ds.^2);
end


        
%% OPTION 2 - If Bayesian inversions are your cup of tea, you'll love this.

if strcmp(invert.inversion_type, 'bayesian') == 1                         %bayesian, whatever your smoothing parameters, VK, laplacain or none  
  
    iterations = invert.iterations;  

    if strcmp(invert.solve_for_fault_size, 'yes') == 1 && strcmp(invert.circular_harmonics, 'yes') == 1
        n_harmonics = 4;
    else
        n_harmonics = 0;
    end

    % Work out how many model parameters you have
    n_model_parameters = total_n_slip_patches;        

         if strcmp(invert.variable_rake, 'yes') ==1                             % If solving for rake
             n_model_parameters = n_model_parameters * 2;
         end

         if strcmp(invert.smoothing, 'none') ~= 1                               % If solving for alpha^2 (one per effective fault strand)
             n_model_parameters = n_model_parameters + n_fault_strands_for_smoothing;
         end

         if strcmp(invert.solve_for_dip, 'yes') == 1                            % If solving for dip (one per effective fault strand)
             n_model_parameters = n_model_parameters + n_fault_strands;
         end

          if strcmp(invert.circular_harmonics, 'yes') == 1                      % If solving for circular harmonics fault size
              n_model_parameters = n_model_parameters + (2 + (n_harmonics*2-1) + 2);  % 2 for xy location of center, then each harmonic has a strength and rotation, then solving for correlation lengths. apart from the first harmonic that doens't have a rotation, because it's a circle
          end
        if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
            n_model_parameters = n_model_parameters + n_InSAR_scenes;           % If solving for InSAR offset 
        end

        if strcmp(invert.solve_for_beta, 'yes') == 1               % If solving for beta (and not setting it at an initial value of 1 and not changing it through the inversion)
            n_model_parameters = n_model_parameters + 1;
        end


    % First, if doing ensemble sampling, set up some things
    if strcmp(invert.ensemble_sampling, 'yes') == 1
        a = 1.1;                                        % This is the scaling distance between which a walker will move between itself and another walker.
        nwalkers = ceil(n_model_parameters/100)*100;    % This rounds the number of walkers up to the nearest 100.   
        iterations = floor(iterations/nwalkers);         % One iteration is defined as when each walker has been perturbed.
        %thinning = 1000;
    else
        nwalkers = 1;
        %thinning = 1;
    end  


    % If solving for fault size 
    if strcmp(invert.solve_for_fault_size, 'yes') == 1
       faultsizecount=1;
       m_on_keep = zeros( n_model_parameters, iterations);
       if strcmp(invert.circular_harmonics, 'yes') == 1                    
           if length(priors.max_circharm_coeffs) ~= n_harmonics
              max_circharm_coeffs  = repmat(priors.max_circharm_coeffs, n_harmonics, 1);
              min_circharm_coeffs  = repmat(priors.min_circharm_coeffs, n_harmonics, 1);
           else
               max_circharm_coeffs  = priors.max_circharm_coeffs;
               min_circharm_coeffs  = priors.min_circharm_coeffs;
           end
           if length(priors.max_circharm_phi) ~= (n_harmonics-1)
              max_circharm_phi = repmat(priors.max_circharm_phi, (n_harmonics-1), 1);
              min_circharm_phi =  repmat(priors.min_circharm_phi, (n_harmonics-1), 1);
           else
               max_circharm_phi = priors.max_circharm_phi;
               min_circharm_phi = priors.min_circharm_phi;
           end
           if length(priors.max_circharm_center) ~= 2   % For x and y
               max_circharm_center = repmat(priors.max_circharm_center, 2, 1);                     % This may seem pointless, but it's set up ready to solve for multiple fault strands
           else
               max_circharm_center = priors.max_circharm_center;
           end
           if length(priors.min_circharm_center) ~= 2
                min_circharm_center = repmat(priors.min_circharm_center, 2, 1);
           else
                min_circharm_center = priors.min_circharm_center; 
           end
           if strcmp(invert.ensemble_sampling, 'yes') == 1  % Initialise walkers
               if strcmp(invert.ensemble_start, 'scatter') == 1 || strcmp(invert.ensemble_start, 'tight') == 1
                   circharm_center_initial(1,:) = fault_length_for_smoothing*rand(1,nwalkers);   % Distribute initial walkers all around the fault plane
                   circharm_center_initial(2,:) =  fault_width_for_smoothing*rand(1,nwalkers);
                    
                   circharm_phi_initial = 2*pi*rand(n_harmonics-1,nwalkers);   % Distribute initial walkers from 0 to 2pi
                   
                   %circharm_coeffs_curr(1,:) = linspace(0, sqrt((fault_length_for_smoothing/2)^2 + (fault_width_for_smoothing/2)^2), nwalkers);
                   %circharm_coeffs_curr(2:n_harmonics,:) = repmat(linspace(0,1000,nwalkers), (n_harmonics-1),1);
                   circharm_coeffs_initial = (max([fault_length_for_smoothing fault_width_for_smoothing])).*rand((n_harmonics),nwalkers);
                   %circharm_coeffs_trial = circharm_coeffs_curr(:,1);
                   %circharm_coeffs_curr(:,1) = [max([fault_length_for_smoothing fault_width_for_smoothing])*2/3; 0; 0; 0]; % Set the first values in circharm_center_curr to such that all the slip patches are on, because this is what we store for the first 1000 iterations
               elseif strcmp(invert.ensemble_start, 'tight') == 1
                   
                   disp('not yet sorted out how to start circharm parameters tight')
                   
               end
           else
               circharm_center_initial = [fault_length_for_smoothing/2; fault_width_for_smoothing/2];       % x,z offset
               circharm_center_step_size_initial = 100*ones(2,1);
               
               circharm_coeffs_initial = [ sqrt((fault_length_for_smoothing/2)^2 + (fault_width_for_smoothing/2)^2); zeros(n_harmonics-1,1)];      % so the circle starts bigger than the fault
               circharm_coeffs_step_size_initial = 100*ones(n_harmonics,1);
               
               circharm_phi_initial = zeros(n_harmonics-1,1);              % Rotation, in radians. 1 rad is ~ 60 degrees.
               circharm_phi_step_size_initial = 0.1*ones(n_harmonics-1,1);  
           end
           
           patchx = 0.5*along_strike_sep_dist+ cumsum([0, repmat(along_strike_sep_dist,1,n_along_strike_patches-1)]);   % Note that solving for fault size ONLY WORKS FOR ONE FAULT currently
           patchx = repmat(patchx,n_down_dip_patches,1);
           patchx = reshape(patchx, 1, []);
           patchz = 0.5*(disloc_model(8,:)+disloc_model(9,:));
           
%             % % Check initial set up
%             [circx,circz]=circharm(circharm_coeffs_initial,circharm_phi_initial,0);     % Set to '1' if want to plot
%             onoffidentifyer = inpolygon(patchx, patchz, (circx+circharm_center_initial(1)), (circz+circharm_center_initial(2)))';
%             mfreeze = zeros((total_n_slip_patches - sum(onoffidentifyer)),1);
%             m_on = [onoffidentifyer; ones((n_model_parameters-total_n_slip_patches), 1)];
%             figure('position', [500, 500, 800, 800]);
%             scatter(patchx(onoffidentifyer==1), patchz(onoffidentifyer==1), 1800, 'gs', 'filled')   % plot x and z of on patches
%             hold on;
%             scatter(patchx(onoffidentifyer==0), patchz(onoffidentifyer==0),1800, 'rs', 'filled')   % plot x and z of off patches
%             plot(circx+circharm_center_initial(1),circz+circharm_center_initial(2));
%             set(gca,'Ydir','reverse')
%             axis equal tight
%             axis([min([circz+circharm_center_initial(1) 0]) max([circz+circharm_center_initial(1)  fault_length_for_smoothing]) min([circx+circharm_center_initial(2) 0]) max([circx+circharm_center_initial(2)  fault_width_for_smoothing])]) 
%             ylabel('Down dip')
%             xlabel('Along strike')
%             title('Slipping area on fault')
       end
    else
       circharm_coeffs_initial = [];
       circharm_phi_initial = [];
       circharm_center_initial = [];
       
       circharm_coeffs_step_size_initial = [];
       circharm_phi_step_size_initial = [];
       circharm_center_step_size_initial = [];
       
       max_circharm_coeffs = [];
       max_circharm_phi = [];
       max_circharm_center = [];
       min_circharm_coeffs = [];
       min_circharm_phi = [];
       min_circharm_center = [];
       
       patchx = [];
       patchz = [];
    end

    if strcmp(invert.ensemble_sampling, 'yes')==1 && strcmp(invert.solve_for_beta, 'yes') == 1  % distribute walkers around starting value
        if strcmp(invert.ensemble_start, 'scatter') == 1
            beta_initial = (10).*rand(n_datasets,nwalkers); % Randomly scatter starting slip parameters within acceptable prior, between 0 and 10 (is 10 reasonable??)
        elseif strcmp(invert.ensemble_start, 'tight') == 1
            beta_initial = linspace(invert.beta_initial/4, invert.beta_initial+invert.beta_initial/4, nwalkers);
        end
    elseif strcmp(invert.ensemble_sampling, 'no')==1 && strcmp(invert.solve_for_beta, 'yes') == 1
        beta_initial = repmat(invert.beta_initial,1,nwalkers);
    elseif strcmp(invert.solve_for_beta, 'no') == 1
        beta_initial = [];
    end

    % Set up initial slip distribution ****************************************
    if n_fault_strands_for_smoothing ~= length(priors.max_slip)
       priors.min_slip = repmat(priors.min_slip, 1, n_fault_strands_for_smoothing);
       priors.max_slip = repmat(priors.max_slip, 1, n_fault_strands_for_smoothing);
    end
    slip_range = priors.max_slip-priors.min_slip;
    max_slip = [];
    min_slip = [];
    for p = 1: n_fault_strands_for_smoothing
        max_slip = [max_slip; priors.max_slip(p)*ones(n_slip_patches_on_each_fault_strand_for_smoothing(p),1)]; % One value per fault patch, in case different fault strands have different priors
        min_slip = [min_slip; priors.min_slip(p)*ones(n_slip_patches_on_each_fault_strand_for_smoothing(p),1)]; % One value per fault patch, in case different fault strands have different priors
    end

    if strcmp(invert.ensemble_sampling, 'yes')==1  % distribute walkers around starting value
        if strcmp(invert.ensemble_start, 'scatter') == 1
            slip_initial = (max_slip-min_slip).*rand(total_n_slip_patches,nwalkers) + min_slip; % Randomly scatter starting slip parameters within acceptable prior
        elseif strcmp(invert.ensemble_start, 'tight') == 1
            walkeroffsets = linspace(-invert.slip_initial/8, invert.slip_initial/8, nwalkers);
            walkeroffsets = repmat(walkeroffsets, total_n_slip_patches,1);
            slip_initial = invert.slip_initial*ones(total_n_slip_patches, nwalkers) + walkeroffsets;
        end
    else
        slip_initial = invert.slip_initial * ones( total_n_slip_patches, nwalkers); 
    end
        if strcmp(invert.ensemble_start, 'tight') ==1 || strcmp(invert.ensemble_sampling, 'no') == 1
            noise_max = 0.01;
            noise_min = -0.01;
            noise = noise_min + (noise_max - noise_min).*rand(total_n_slip_patches,nwalkers);
            slip_initial = slip_initial + noise;
                for j = 1: total_n_slip_patches
                    if slip_initial(j) > max_slip(j)
                        slip_initial(j) = 2*max_slip(j) - slip_initial(j);
                    end
                    if slip_initial(j) < min_slip(j)
                        slip_initial(j) = (2* min_slip(j)) - slip_initial(j);     
                    end
                end  
        end

         slip_initial(onoffidentifyer==0) = 0;
        
    disp('Initial slip distribution set up. Onwards!')
    disp('     ')


    % Housekeeping.............................................................
    store_number= 0;
    %total_rejections = zeros(1,invert.iterations);
    rejection_target = 0.766;  %  Roberts et al 1997
    m_keep = zeros(n_model_parameters, iterations, nwalkers);
    M0_keep = zeros(1, ceil(invert.iterations));
    M0_likelihood_keep = zeros(1,ceil(invert.iterations));
    logprior_keep = zeros(n_fault_strands, ceil(invert.iterations));
    logposterior_keep = zeros(1, ceil(invert.iterations));
    onoffidentifyerkeep = zeros(total_n_slip_patches, 1);
    det_sigma_s_keep = zeros(n_fault_strands_for_smoothing,invert.iterations);
    logL_keep = zeros(1, ceil(invert.iterations));
    transdimensional_prior_ratio = 1;

    det_sigma_s = [];           % Not necessary if not doing a VK inversion
    det_sigma_s_ratio = 1;      % This will only be updated if solving for faultsize i.e. if changing the sigma_s between iterations
    sigma_s_master = [];        % Not necessary if not doing a VK inversion
    inv_sigma_s_master = [];    % Not necessary if not doing a VK inversion
    L =[];                      % Not necessary if not doing a laplacian inversion

        if strcmp(priors.slip_prior, 'gaussian') + strcmp(priors.slip_prior, 'logarithmic')  > 0
            b_curr = zeros(total_n_slip_patches, 1);
            max_b_trial = log10(max_slip);
            min_b_trial = log10(min_slip);
        end

    offset_step_size_initial = [];
    max_offset = [];
    min_offset = [];
       
    if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1 
        if n_InSAR_scenes ~= length(priors.min_offset)
            priors.min_offset = repmat(priors.min_offset, n_InSAR_scenes, 1);
            priors.max_offset = repmat(priors.max_offset, n_InSAR_scenes, 1);
        end
        min_offset = priors.min_offset;
        max_offset = priors.max_offset;
        if strcmp(invert.ensemble_sampling, 'yes')==1
            if strcmp(invert.ensemble_start, 'scatter') == 1
                for i = 1:n_InSAR_scenes
                    offset_initial(InSAR_identifyer==i,:) = repmat(((priors.max_offset-priors.min_offset).*rand(1,nwalkers) + priors.min_offset),sum(InSAR_identifyer==i),1); % Randomly scatter starting slip parameters within acceptable prior
                end
            elseif strcmp(invert.ensemble_start, 'tight') == 1  % Distribute walkers around starting value
                for i = 1:n_InSAR_scenes
                    walkeroffsets = linspace(-0.1, 0.1, nwalkers);
                    walkeroffsets = repmat(walkeroffsets, sum(d_identifyer==1),1);
                    offset_initial = invert.offset_initial+walkeroffsets;
                end
                offset_initial = [offset_initial; zeros(sum(d_identifyer==2), nwalkers)];   % zeros on bottom for GPS
            end
        elseif strcmp(invert.ensemble_sampling, 'no')==1
            offset_step_size_initial = 0.01*ones(n_InSAR_scenes,1);    %meters
            offset_modelparameter_initial = invert.offset_initial*ones(n_InSAR_scenes,1);
            offset_initial(d_identifyer==1,:) = invert.offset_initial*ones(sum(d_identifyer==1), nwalkers);       % one row per InSAR dataset, one column per walker number
            offset_initial = [offset_initial; zeros(sum(d_identifyer==2), nwalkers)];   % zeros on bottom for GPS
        end
    else
        offset_modelparameter_initial = []; 
        offset_initial = []; 
    end

    % Sort out G matrix for solving for dip.................................... 
    
            if strcmp(invert.solve_for_dip, 'yes') == 1
                if n_fault_strands ~= length(priors.min_dip)        
                    priors.min_dip = repmat(priors.min_dip, 1, n_fault_strands_for_smoothing);
                    priors.max_dip = repmat(priors.max_dip, 1, n_fault_strands_for_smoothing);
                end
                G_ss_curr = zeros(n_data, total_n_slip_patches);
                G_ds_curr = zeros(n_data, total_n_slip_patches);    
               % dip_initial = ((max_dip - min_dip).*rand(n_fault_strands, 1) + min_dip);  % one per fault strand (NOT FAULT STRAND PER SMOOTHING)
                dip_initial = invert.dip_initial;
                max_dip = priors.max_dip;
                min_dip = priors.min_dip;
                    for k = 1 : n_fault_strands       % select right G matrix for that dip. recall each fault strand can have a different dip
                        G_ss_curr(:, first_patch_in_strand(k):last_patch_in_strand(k)) = interp3(1:n_slip_patches_on_each_fault_strand(k), 1:n_data, dip_LUT, G_ss(:,first_patch_in_strand(k):last_patch_in_strand(k),:), 1:n_slip_patches_on_each_fault_strand(k), 1:n_data, dip_initial(k));
                        G_ds_curr(:, first_patch_in_strand(k):last_patch_in_strand(k)) = interp3(1:n_slip_patches_on_each_fault_strand(k), 1:n_data, dip_LUT, G_ds(:,first_patch_in_strand(k):last_patch_in_strand(k),:), 1:n_slip_patches_on_each_fault_strand(k), 1:n_data, dip_initial(k));
                    end 
            elseif strcmp(invert.solve_for_dip, 'no') == 1
                dip_initial = [];
                G_ss_curr = G_ss;
                G_ds_curr = G_ds;
                G_ss_temp = G_ss;
                G_ds_temp = G_ds;
                max_dip = [];
                min_dip = [];
            end


    % Variable rake........................................ ...................     
    max_rake = [];
    min_rake = [];
        
          if strcmp(invert.variable_rake, 'yes') == 1

              if n_fault_strands_for_smoothing ~= length(priors.min_rake)
                  priors.min_rake = repmat(priors.min_rake, 1, n_fault_strands_for_smoothing);
                  priors.max_rake = repmat(priors.max_rake, 1, n_fault_strands_for_smoothing);
              end
              rake_range = priors.max_rake-priors.min_rake;
              
              for p = 1: n_fault_strands_for_smoothing
                  max_rake = [max_rake; priors.max_rake(p)*ones(n_slip_patches_on_each_fault_strand_for_smoothing(p),1)]; % One value per fault patch, in case different fault strands have different priors
                  min_rake = [min_rake; priors.min_rake(p)*ones(n_slip_patches_on_each_fault_strand_for_smoothing(p),1)]; % One value per fault patch, in case different fault strands have different priors
              end

               rake_initial = disloc_model(5, :)';

               if strcmp(invert.ensemble_sampling, 'yes') == 1  % distribute walkers around starting value
                    if strcmp(invert.ensemble_start, 'scatter') == 1
                        for i = 1:n_fault_strands_for_smoothing
                            rake_initial(first_patch_in_strand_for_smoothing(i):last_patch_in_strand_for_smoothing(i), 1:nwalkers) = (priors.max_rake(i)-priors.min_rake(i)).*rand(n_slip_patches_on_each_fault_strand_for_smoothing(i),nwalkers) + priors.min_rake(i);
                        end
                    elseif strcmp(invert.ensemble_start, 'tight') == 1
                        for i = 1:n_fault_strands_for_smoothing
                            walkeroffsets = linspace(-10,10,nwalkers);
                            walkeroffsets = repmat( walkeroffsets, n_slip_patches_on_each_fault_strand_for_smoothing(i), 1);
                            rake_initial(first_patch_in_strand_for_smoothing(i):last_patch_in_strand_for_smoothing(i), 1:nwalkers) = mean([priors.min_rake priors.max_rake]) + walkeroffsets;
                            noise_max = 3;
                            noise_min = -3;
                            noise = noise_min + (noise_max - noise_min).*rand(total_n_slip_patches,nwalkers);
                            rake_initial = rake_initial + noise;
                        end
                    end
               end

               G_curr = G_ss_curr*diag(cosd(rake_initial)) +  G_ds_curr*diag(sind(rake_initial));  
           
         else
             rake_initial = [];
             G_curr = G;  % Since not solving for variable rake, G_curr = G_temp = G, for the whole iteration, which is just the G calculated using the rake given in fault.fault_descriptor_file
             G_temp = G;
          end

         G_sensitivity_test = G_curr;   % this gets updated later. but we want to start it the right size, as we only update the rows for the slip patches that are on




    % Step sizes and such......................................................

    %patch_step_sizes = (invert.step_size / n_model_parameters) * ones(total_n_slip_patches, 1) * 100;  
    patch_step_sizes_initial = invert.step_size* ones(total_n_slip_patches, 1);
    rake_step_size_initial = [];
    if strcmp(invert.variable_rake, 'yes') ==1
        rake_step_size_initial = 0.5*ones(total_n_slip_patches, 1);
    end
    if strcmp(invert.solve_for_dip, 'yes') == 1
        dip_step_size_initial = ones(n_fault_strands, 1);
    else
        dip_step_size_initial = [];
    end
    %probability_target = invert.probability_target_initial / n_model_parameters * 100;
    probability_target = invert.probability_target_initial;
    sens_test_number = 1;
    sens_test = [100:100:900, 1000:500:5000, 6000:1000:10000, 20000:10000:invert.iterations];
    store_number_at_sens_test = zeros(1, length(sens_test));    % preallocate for speed    
    sens_test(end+1) = 0;      % add a 0 on the end, because there's a counter to let you know which sensitivity test you're on, and so the last sensitivity test is one bigger than the length of sens_test, and matlab doesn't like checking an empty matrix entry. so pad it out.
    %step_sizes_keep = zeros(n_, length(sens_test));      
    SEM = zeros(total_n_slip_patches, (length(sens_test)-1));
    CI_low = zeros(total_n_slip_patches, (length(sens_test)-1));
    CI_high = zeros(total_n_slip_patches, (length(sens_test)-1));
    %rej_ratio_since_last_sens_test = 0;     
    sens_test_interval = diff([0 sens_test]);     
    %sens_test_interval(end) = []; 
    step_sizes_keep = zeros(n_model_parameters, length(sens_test));
    rej_ratio_save = zeros(1,length(sens_test));


    % *********** % ************ % ********** % ********** % ********** % *****
    %  Unsmoothed bayesian specific housekeeping...............................
    % *********** % ************ % ********** % ********** % ********** % *****   

    if strcmp(invert.inversion_type, 'bayesian') == 1 && strcmp(invert.smoothing, 'none') == 1 
            if strcmp(invert.pad_edges_with_zeros, 'yes') == 1
               slip_initial(1: n_down_dip_patches) = 0;                             % first column = 0
               slip_initial((end-n_down_dip_patches): end) = 0;                     % last column = 0
               slip_initial(n_down_dip_patches: n_down_dip_patches: end) = 0;       % bottom row = 0
            end
        logL_curr = calc_loglikely(slip_initial, d, G_curr, inv_sigma_d, offset_curr, beta_initial);   %  not exp((-1/2) * ( (d - d_hat).' * inv_sigma_d * (d - d_hat) )) because we're taking logs; 
        logposterior_curr = logL_curr;
        inv_sigma_s = zeros(total_n_slip_patches);       % not necessary if not doing a VK inversion
        alpha2_curr = [];        % not necessary if not doing a VK inversion
        alpha2_step_size_initial = [];   % not necessary if not doing a VK inversion

        ratio_keep = zeros(ceil(invert.iterations),1);
        H= [];
        predominant_faulting_style = [];
    end



    % *********** % ************ % ********** % ********** % ********** % *****
    %  Smoothing bayesian specific housekeeping...............................
    % *********** % ************ % ********** % ********** % ********** % *****   
    if strcmp(invert.smoothing, 'VK') || strcmp(invert.smoothing, 'laplacian') || strcmp(invert.smoothing, 'minimumnorm') == 1

        % Pre-locate some memory...........................................
        pre_prior_rejections = zeros(1,invert.iterations);
        prior_rejections = zeros(1,invert.iterations);
        likelihood_rejections = zeros(1,invert.iterations);
        
        % Initialise prior hyperperameter, alpha^2
        alpha2_keep = zeros(n_fault_strands_for_smoothing,ceil(invert.iterations));
        alpha2_modelparameter_initial = invert.alpha2_initial;
            if n_fault_strands_for_smoothing ~= length(priors.min_alpha2)
                priors.min_alpha2 = repmat(priors.min_alpha2, 1, n_fault_strands_for_smoothing);
                priors.max_alpha2 = repmat(priors.max_alpha2, 1, n_fault_strands_for_smoothing);
                invert.alpha2_initial = repmat(invert.alpha2_initial, 1, n_fault_strands_for_smoothing);
                alpha2_modelparameter_initial = invert.alpha2_initial';
            end
            if strcmp(priors.alpha2_prior, 'logarithmic') == 1
                %alpha2_modelparameter_initial = repmat(log10(invert.alpha2_initial)', n_fault_strands_for_smoothing, nwalkers);
                %alpha2_modelparameter_initial = invert.alpha2_initial;
                if priors.min_alpha2 == 0
                   min_alpha2 = 0.001 * ones(length(priors.min_alpha2), 1); 
                else
                   min_alpha2 = priors.min_alpha2;
                end
                max_alpha2modelparameter = log10(priors.max_alpha2)';
                min_alpha2modelparameter = log10(min_alpha2)';
            else
                alpha2_modelparameter_initial = repmat(invert.alpha2_initial, n_fault_strands_for_smoothing,nwalkers);
                max_alpha2modelparameter = priors.max_alpha2';
                min_alpha2modelparameter = priors.min_alpha2';
            end

       if strcmp(invert.ensemble_sampling, 'yes') == 1  % distribute walkers around starting value
            if strcmp(invert.ensemble_start, 'scatter') == 1
                alpha2_modelparameter_initial = (max_alpha2modelparameter-min_alpha2modelparameter).*rand(n_fault_strands_for_smoothing,nwalkers)+min_alpha2modelparameter;
            elseif strcmp(invert.ensemble_start, 'tight') == 1
                for p = 1: n_fault_strands_for_smoothing
                    alpha2_modelparameter_initial(p,:) = invert.alpha2_initial(p)+ linspace( -invert.alpha2_initial(p)/2,invert.alpha2_initial(p)/2, nwalkers);
                    if strcmp(priors.alpha2_prior, 'logarithmic') == 1
                        alpha2_modelparameter_initial = log10(alpha2_modelparameter_initial);
                    end
                end
            end
       end
       
       if size(alpha2_modelparameter_initial,2) > size(alpha2_modelparameter_initial,1)
           alpha2_modelparameter_initial = alpha2_modelparameter_initial';
       end

       if strcmp(invert.ensemble_sampling, 'yes') == 1 && n_fault_strands_for_smoothing>1
           disp('HEY - check how alpha^2 is initialised for more than one fault strand with ensemble sampling')
       end

       if strcmp(priors.alpha2_prior, 'logarithmic') == 1
            alpha2_initial = 10.^alpha2_modelparameter_initial;
        else
            alpha2_initial = alpha2_modelparameter_initial;
        end

        if strcmp(invert.smoothing, 'minimumnorm') == 1 && any(alpha2_initial<1) && strcmp(priors.alpha2_prior, 'logarithmic') == 0
            alpha2_initial = ones(n_fault_strands_for_smoothing,1);
            disp('alpha2_initial was less than 1 so I''ve changed it to 1');
        end
        alpha2_step_size_initial = (invert.alpha2_step_size.*ones(1,n_fault_strands_for_smoothing))'; %repmat(invert.alpha2_step_size, n_fault_strands_for_smoothing,1);
        if strcmp(invert.smoothing, 'minimumnorm') == 1 && any(alpha2_step_size_initial<0.01) ==1 && strcmp(priors.alpha2_prior, 'logarithmic') == 0
            alpha2_step_size_initial = repmat(0.1, n_fault_strands_for_smoothing,1);
            disp('alpha2_step_size_initial was too small so I changed it to 0.01');
        end
        %prior_ratio_keep = zeros(ceil(invert.iterations), 1);
        alpha2_curr = alpha2_initial;
        

        %.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.
        % Laplacian Bayesian specific housekeeping %.%.%.%.%.%.%.%.%.%.%.%.%.%.
        %.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.

                    a_as_modelparameter_initial = [];
                    a_dd_modelparameter_initial = [];
                    a_as_step_size = [];
                    a_dd_step_size = [];
                    max_a_as_modelparameter = [];
                    min_a_as_modelparameter = [];
                    max_a_dd_modelparameter = [];
                    min_a_dd_modelparameter = [];
        
        if strcmp(invert.smoothing, 'laplacian') == 1 
                r = down_dip_sep_dist./along_strike_sep_dist;  %   r is ratio of patch length in 1st dim over length in 2nd dim
                free_edge= 1;   % 1 = top free
                
                for i = 1: n_fault_strands_for_smoothing   % If we're not smoothing across the faults, then we'll treat them all separately, so we'll store them all in a big diagonal matrix.                
                    L_temp = laplace_fault(n_down_dip_patches_for_smoothing(i),n_along_strike_patches_for_smoothing(i),r(i),free_edge);  % thanks to andy for this one
                    L = blkdiag(L,L_temp);
                    clear L_temp
                end
                
               % THESIS CORRECTIONS - trying out all free edges
               L =[];
                for i = 1: n_fault_strands_for_smoothing   % If we're not smoothing across the faults, then we'll treat them all separately, so we'll store them all in a big diagonal matrix.                
                    L_temp = laplace_fault_allfree(n_down_dip_patches_for_smoothing(i),n_along_strike_patches_for_smoothing(i),r(i),free_edge);  % thanks to andy for this one
                    L = blkdiag(L,L_temp);
                    clear L_temp
                end 
                free_edge = 'all_free';

                logprior_curr = calc_logprior_laplacian( slip_initial(:,1), L, alpha2_initial(:,1), n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master);
                predominant_faulting_style = [];
                H= [];                      % Not necessary if not doing VK inversion
                inv_sigma_s = [];           % Not necessary if not doing VK inversion
                det_sigma_s = [];           % Not necessary if not doing VK inversion

        %.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.
        % Minimum norm Bayesian specific housekeeping %.%.%.%.%.%.%.%.%.%.%.%.%.
        %.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.

        elseif strcmp(invert.smoothing, 'minimumnorm') == 1
                for q = 1:n_fault_strands_for_smoothing
                    s_curr = slip_initial(first_patch_in_strand_for_smoothing(q):last_patch_in_strand_for_smoothing(q),1);
                    logprior_curr(q,1) =  ( -1/(2*alpha2_initial(q,1).^2) * sum((s_curr'*s_curr)));
                end
                H= [];

        %.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%    
        % VK Bayesian specific housekeeping %.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%
        %.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%.%
          
        
        elseif strcmp(invert.smoothing, 'VK') == 1

                % The correlation lengths depend on the predominant faulting style - which may not be the same for all fault strands.
                if length(priors.predominant_faulting_style) ~= n_fault_strands_for_smoothing
                    predominant_faulting_style = repmat(priors.predominant_faulting_style,[1,n_fault_strands_for_smoothing]);   
                else
                    predominant_faulting_style = priors.predominant_faulting_style;
                end

                % Calculate correlation lengths
                for n = 1 : n_fault_strands_for_smoothing 
                    if strcmp(predominant_faulting_style(n), 'ss') == 1
                        a_as(n) =  1860 + 0.34*(fault_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
                        a_dd(n) =  -390 + 0.44 * (fault_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002) 
                    elseif strcmp(predominant_faulting_style(n), 'ds') == 1 
                        a_as(n) =  1100 + 0.31*(fault_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
                        a_dd(n) =  580 + 0.35* (fault_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002) 
                    end
                end
               
                if strcmp(invert.solve_for_fault_size, 'yes') == 1         % If solving for fault size we'll solving for a_as and a_dd, drawing from their distributions
                    a_as_modelparameter_initial = 0;
                    a_dd_modelparameter_initial = 0;
                    a_as_step_size = 0.05;
                    a_dd_step_size = 0.05;
                    max_a_as_modelparameter = 1;           % To draw from a normal distribution with a random walk we use the error function trick (can't do random walk with randn). The erfinv function only exists between -1 and 1
                    min_a_as_modelparameter = -1;
                    max_a_dd_modelparameter = 1;
                    min_a_dd_modelparameter = -1;
                    a_as_keep = zeros(invert.iterations,1);
                    a_dd_keep = zeros(invert.iterations,1);
                    if strcmp(predominant_faulting_style, 'ss') == 1
                        % Strike slip parameters
                        a_as_ss_mean = 1860 + 0.34*fault_length_for_smoothing;
                        a_as_ss_std = sqrt(1120^2+0.03^2*fault_length_for_smoothing^2);          % Errors from (Mai and Beroza, 2002), Table 2, Pg 14
                        a_dd_ss_mean = -390 + 0.44 *fault_width_for_smoothing;                   % Propagation of errors from talbe on Wikipedia (shh)  https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                        a_dd_ss_std = sqrt(470^2+0.04^2*fault_width_for_smoothing^2);
                    elseif strcmp(predominant_faulting_style, 'ds') == 1
                        % Dip slip parameters
                        a_as_ds_mean =  1100 + 0.31*fault_length_for_smoothing;
                        a_as_ds_std = sqrt(400^2+0.01^2*fault_length_for_smoothing^2);
                        a_dd_ds_mean =  580 + 0.35*fault_width_for_smoothing;
                        a_dd_ds_std = sqrt(240^2+0.01^2*fault_width_for_smoothing^2);
                    end
                else
                    %a_as_modelparameter_initial = [];
                    %a_dd_modelparameter_initial = [];
                    %a_as_step_size = [];
                    %a_dd_step_size = [];
                    %max_a_as_modelparameter = [];
                    %min_a_as_modelparameter = [];
                    %max_a_dd_modelparameter = [];
                    %min_a_dd_modelparameter = [];
                end
                %proposal_ratio = 1;     % this is the ratio that changes is the dimension of the problem changes (i.e. if a transdimensional step is taken).

                % Calculate scaled distances between all fault patches
                [r_over_a, angle_between_patches] = calc_scaled_dist( n_fault_strands_for_smoothing, disloc_model, a_as, a_dd, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master, along_strike_sep_dist, n_along_strike_patches, n_down_dip_patches, fault_strand_togetherness);
                H_as = 0.71;                                  % Along-strike direction-dependent Hurst number  (Mai and Beroza, 2002) 
                H_dd = 0.77;                                  % Down-dip direction-dependent Hurst number  (Mai and Beroza, 2002) 
      
                % Calculate direction dependent H
                H = [];
                for n = 1: n_fault_strands_for_smoothing
                    H_temp = calc_direction_dependent_H( H_dd, H_as, angle_between_patches(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n))); % along strike is purely H_as, down dip is purely H_dd, and everything else is somewhere inbetween, scaled geometrically
                    H = blkdiag(H,H_temp);
                    H_temp = [];
                end    

                % Calculate sigma_s
                sigma_s = [];
                for i = 1: n_fault_strands_for_smoothing
                    sigma_s_temp = calc_sigma_s( r_over_a(first_patch_in_strand_for_smoothing_master(i):last_patch_in_strand_for_smoothing_master(i),first_patch_in_strand_for_smoothing_master(i):last_patch_in_strand_for_smoothing_master(i)), H(first_patch_in_strand_for_smoothing_master(i):last_patch_in_strand_for_smoothing_master(i),first_patch_in_strand_for_smoothing_master(i):last_patch_in_strand_for_smoothing_master(i)));       % NOTE that this is actually a SEPARATE SIGMA S MATRIX FOR EACH FAULT STRAND, just stored in one big matrix
                    sigma_s = blkdiag(sigma_s,sigma_s_temp);
                    clear sigma_s_temp
                end
                if strcmp(invert.add_correlation_matrix_stabiliser, 'yes') == 1
                    sigma_s = sigma_s + diag( 0.01*ones(total_n_slip_patches,1));
                end
                sigma_s_master = sigma_s;
                det_sigma_s = [];
                inv_sigma_s = [];

                % Calculate det_sigma_s and inv_sigma_s
                    for n = 1: n_fault_strands_for_smoothing            % NOTE that this is actually a SEPARATE invSIGMA S MATRIX FOR EACH FAULT PATCH, just stored in one big matrix
                        det_sigma_s(i,1) = det(sigma_s(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n), first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)));       
                        inv_sigma_s_temp(:, :) = inv(sigma_s(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n), first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)));
                        inv_sigma_s = blkdiag(inv_sigma_s, inv_sigma_s_temp);
                        inv_sigma_s_temp = [];
                    end
                    clear inv_sigma_s_temp
                    det_sigma_s_keep(:,1) = det_sigma_s;
                    inv_sigma_s_master = inv_sigma_s;
                    
                    if det_sigma_s == 0
                       error('det_sigma_s is 0, which usually means you have too many slip patches, which usually means the inversion does not work'); 
                    end
                    
        end
        
        % Moment regularisation ...............................................
        if strcmp(invert.regularise_moment, 'yes') == 1
            
            M0_curr = sum((elastic_params.mu_okada) * spatial_model2column .* spatial_model3column .* (slip_initial(:,1)));
            M0_likelihood_curr = normpdf( M0_curr, data.seismic_moment, data.moment_std);
            
            while M0_likelihood_curr == 0
                prompt = 'M0 likelihood is 0. Pick a new value of slip_initial, type it here and press enter: ';
                disp('Thanks!');
                new_slip_initial = input(prompt);
                slip_initial = new_slip_initial + noise;
                M0_curr = sum((elastic_params.mu_okada) * spatial_model2column .* spatial_model3column .* (slip_initial));
                M0_likelihood_curr = normpdf( M0_curr, data.seismic_moment, data.moment_std);
                
                logposterior_curr = logposterior_curr + log(M0_likelihood_curr);
                
            end
            
        elseif strcmp(invert.regularise_moment, 'no') == 1
            M0_likelihood_curr = [];
        end
        
    end
    
    

    % Set up m matrices
    % 1 = slip, 2 = alpha^2, 3 = rake, 4 = dip, 5 = offset, 6 = beta, 7 = circhram coeffs, 8 = circharm phi, 9 = circharm centre, 10 = a_as model parameter, 11 = a_dd model parameter
    m_initial = [slip_initial; alpha2_modelparameter_initial; rake_initial; dip_initial; offset_modelparameter_initial; beta_initial; circharm_coeffs_initial; circharm_phi_initial; circharm_center_initial; a_as_modelparameter_initial; a_dd_modelparameter_initial];    % Rows = model parameters, Columns = walkers
    m_curr = m_initial;
    m_trial = m_curr(:,1);
    m_identifyer_master = [ 1* ones(size(slip_initial,1),1); 2*ones(size(alpha2_modelparameter_initial,1),1); 3*ones(size(rake_initial,1),1); 4*ones(size(dip_initial,1),1); 5*ones(size(offset_modelparameter_initial,1),1); 6*ones(size(beta_initial,1),1); 7*ones(size(circharm_coeffs_initial,1),1); 8*ones(size(circharm_phi_initial,1),1); 9*ones(size(circharm_center_initial,1),1); 10*ones(size(a_as_modelparameter_initial,1),1); 11*ones(size(a_dd_modelparameter_initial,1),1)];
    m_max = [max_slip; max_alpha2modelparameter; max_rake; max_dip; max_offset; max_beta; max_circharm_coeffs; max_circharm_phi; max_circharm_center; max_a_as_modelparameter; max_a_dd_modelparameter];
    m_min = [min_slip; min_alpha2modelparameter'; min_rake; min_dip; min_offset; min_beta; min_circharm_coeffs; min_circharm_phi; min_circharm_center; min_a_as_modelparameter; min_a_dd_modelparameter];
    if strcmp(invert.ensemble_sampling, 'no') == 1
        step_sizes = [patch_step_sizes_initial; alpha2_step_size_initial; rake_step_size_initial; dip_step_size_initial; offset_step_size_initial; beta_step_size_initial; circharm_coeffs_step_size_initial;  circharm_phi_step_size_initial; circharm_center_step_size_initial; a_as_step_size; a_dd_step_size];
        %step_sizes_master = step_sizes;
    end
    m_identifyer = m_identifyer_master;     % If not solving for fault size, they will always be equal to each other. But if we remove patches, this will be updated
    circharmparameters = zeros(n_model_parameters,1);                       % Same length as the number of model parameters that are on
    circharmparameters(m_identifyer==7|m_identifyer==8|m_identifyer==9|m_identifyer==10|m_identifyer==11) = 1;
    m_on = [onoffidentifyer; ones(n_fault_strands_for_smoothing,1); onoffidentifyer; ones(size(dip_initial,1),1); ones(size(offset_modelparameter_initial,1),1); ones(size(beta_initial,1),1); ones(sum(circharmparameters), 1)];
    m_on_keep(:, 1) = m_on;
     mcircharmfreeze = m_initial(circharmparameters==1);
    
    % Calculate probabilities for initial slip distribution - this is so your first trial has a previous 'curr' to compare to 
    singularflag = 0;
    for ii = 1:nwalkers
        if strcmp(invert.smoothing, 'VK') == 1
           [logprior_curr(:,ii), singularflag] = calc_logprior_VK(m_initial(m_identifyer==1,ii), inv_sigma_s, alpha2_initial(:,ii), n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master);        % NOTE: simple means we don't calculate anything other than the exponent, because for non varying H and a, this is the same every time.
        end
        logL_curr(1,ii) = calc_loglikely(m_initial(m_identifyer==1,ii), d, G_curr, inv_sigma_d, offset_initial, m_initial(m_identifyer==6,ii));
        logposterior_curr(1,ii) =  log(sum( (2*pi*alpha2_initial(:,ii)).^(-n_slip_patches_on_each_fault_strand_for_smoothing/2))) + sum(logprior_curr(1,ii)) + logL_curr(:,ii);                      % if we hadn't taken logs, we'd multiply them. but we calculated the log probability, so we add them. log(AB) = log(A) + log(B).
    end
    
    if strcmp(invert.solve_for_beta, 'yes') == 1
       logposterior_curr = logposterior_curr * (length(d)/2)*log(2*pi*m_initial(m_identifyer==6,ii)^2);
    end
    
    if strcmp(invert.regularise_moment, 'yes') == 1
        logposterior_curr = logposterior_curr + log(M0_likelihood_curr);        % VK logposterior_curr
    end
                 
    disp('.') 

    % Load old MCMC chain, if you want to
    if strcmp(invert.load_old_MCMC_chain, 'no') == 0  % If there is a name of a file given in 'load_old_MCMC_chain' then we'll load this and use this

           disp('Loading old MCMC chain...')
           load(invert.load_old_MCMC_chain);
           
           % Preallocate some space so the inversion takes less time
           store_number = size(slip_keep,2);    % Need to add to the end of the matrixes. but some may have been removed from the burn-in
           rake_keep(:,(store_number+1:store_number+iterations)) = 0;
           logL_keep(:,(store_number+1:store_number+iterations)) = 0;
           M0_keep(1,(store_number+1:store_number+iterations)) = 0;
           M0_likelihood_keep(1,(store_number+1:store_number+iterations)) = 0;
           slip_keep(:, (store_number+1:store_number+iterations)) = 0;             % One column for each iteration, containing slip values for each patch (number of rows = number of slip patches)
           logposterior_keep(1,(store_number+1:store_number+iterations)) = 0;  % One column for each iteration, containing posterior for each patch (number of rows = number of slip patches) THIS ISN'T RIGHT FOR MULTIPLE FAULT STRANDS so I've commented it out .
           if strcmp(invert.smoothing, 'none') == 1
               ratio_keep((store_number+1:store_number+iterations),1) = 0;
           end
           if strcmp(invert.smoothing, 'VK') == 1
               alpha2_keep(:,(store_number+1:store_number+iterations)) = 0;     % One column for each iteration, containing alpha2_squared value for each iteration
%                if strcmp(invert.solve_for_correlation_length, 'yes') == 1
%                    a_as_keep(1:n_fault_strands_for_smoothing,(store_number+1:store_number+iterations)) = 0;
%                    a_dd_keep(1:n_fault_strands_for_smoothing,(store_number+1:store_number+iterations)) = 0;
%                end
           end
           if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
               offset_keep(:,(store_number+1:store_number+iterations)) = 0;
           end
           
           housekeeping.save_name = ['loadoldchain', housekeeping.save_name];
    end


    %% Calculate resolution matrices


    if strcmp(display.plot_resolution_matrix, 'yes') == 1
        GT = G';
        Gg = inv(GT * inv_sigma_d * G + inv_sigma_s) * GT * inv_sigma_d;
        resolution_matrices = Gg * G;
        figure
        imagesc(resolution_matrices)
        colorbar
        caxis([0 1])
        title('Resolution matrix')

        patch_resolution = diag(resolution_matrices);

        faults = disloc_model;
        faults(6,:) = patch_resolution;
        figure
        doplot3d(faults', 'jet')
        title('Resolution')
        colorbar
    end


    %%  ************************************************************************ 
    % Simulated annealing *****************************************************
    % ************************************************************************* 

    if strcmp(invert.simulated_annealing_start, 'yes') == 1
        
        disp('Performing simulated annealing to find best start values.')
        simulated_annealing_addin.m

    end
        
    
    
%% ************************************************************************ 
% Quick check! A stitch in time saves three and a half years of misery ****
% *************************************************************************    
    
if strcmp(invert.quickcheck, 'yes') == 1
    check_setup
    disp('   ');
    disp('Is all your data and fault in the right place?');
    disp('If so, click anywhere on the figure to continue.');
    disp('If not you''d best sort that out...');
    disp('   ');
    waitforbuttonpress;
end

    %% ************************************************************************ 
    % LOOP LOOP LOOP **********************************************************
    % ************************************************************************* 

    burn_in_remove_number = 1;
    true_trials = 0;

    if strcmp(invert.smoothing, 'none')
        disp('Starting MCMC - Unsmoothed Bayesian inversion') 
    elseif strcmp(invert.smoothing, 'VK')
        disp('Starting MCMC - Von Karman regularised Bayesian inversion') 
    elseif strcmp(invert.smoothing, 'laplacian')
        disp('Starting MCMC - Laplacian smoothed Bayesian inversion') 
    end 

    iifaultsize = 1;    % The first values in circharm_center_curr and circharm_coeffs_curr are set up so that all patches on the fault is on. This is what will be saved for the first 1000 iterations because the slipping area is changed.
    if strcmp(invert.solve_for_fault_size, 'yes') ==1
        faultsizechangefreq = 2;
        faultsizechangestart = 2;
    else
        faultsizechangefreq = [];
        faultsizechangestart = invert.iterations;
    end
    walkersubsets = [1 ceil( nwalkers/2)];                                                                              % On row 1,  Subset 1 = 1: half the number of walkers                
    walkersubsets = [walkersubsets; max(walkersubsets)+1 nwalkers];                                                     % On row 2,  Subset 2 = half the number of walkers + 1 : nwalkers
    walkerothersubsetidentifyer = [2*ones( walkersubsets(1,2),1);  ones(walkersubsets(2,2)-walkersubsets(2,1)+1,1)];    % We want to perform stretch/walk move with a walker from the other subset. This matrix tells us the 'other' subset for each walker.
    n_slip_patches_ON_on_each_fault_strand_for_smoothing_save = zeros(iterations, n_fault_strands_for_smoothing);
    
    for i = 1:iterations

        for ii = 1: nwalkers

            store_number = store_number+1;  % this is different from i, in the case of ensemble sampling
            
            % Either perturb the paramterisation . . .            
            if strcmp(invert.solve_for_fault_size, 'yes') == 1 && store_number >= faultsizechangestart && rem(store_number,faultsizechangefreq) == 0 || strcmp(invert.solve_for_fault_size, 'yes') == 1 && faultsizechangestart == 1 && store_number == 1   % Update which patches are on and off, at chosen frequency, once left burn in
                iifaultsize = iifaultsize+1;
                faultsizecount = faultsizecount+1;
                if strcmp(invert.ensemble_sampling, 'yes') == 1
                    if strcmp(invert.ensemble_move_style, 'stretch') == 1
                        [ m_trial(circharmparameters==1) ] = do_ensemble_stretch_move( m_curr(circharmparameters==1,:), ii, otherwalker, z);
                    elseif strcmp(invert.ensemble_move_style, 'walk') == 1
                        [ m_trial(circharmparameters==1) ] = do_ensemble_walk_move( m_curr(circharmparameters==1), ii, otherwalkers);
                    end    
                elseif strcmp(invert.ensemble_sampling, 'no') == 1
                    m_curr(circharmparameters==1) = mcircharmfreeze;                % if everything's been rejected, revert back to the patches that were on before the most recent change of slipping area. if not this has been updated to the current circharmparameters, after second metropolis step
                    m_trial(circharmparameters==1) = m_curr(circharmparameters==1) + ((step_sizes(circharmparameters==1) - (-step_sizes(circharmparameters==1))) .*rand(sum(circharmparameters==1),1) - step_sizes(circharmparameters==1));
                    m_curr(circharmparameters==1) = m_trial(circharmparameters==1); % update m_curr for these 500 iterations, even if trial i=500 is rejected
                end
                
            else            % . . . Or perturb the parameters 

                if strcmp(invert.ensemble_sampling, 'yes') == 1
                    % Pick another walker with which to update current walker
                    othersubset = walkerothersubsetidentifyer(ii);
                    if strcmp(invert.ensemble_move_style, 'stretch') == 1                                                                       % Randomly select which other walker will be participating in this move
                        %otherwalker = ceil(walkersubsets(othersubset,1) + (walkersubsets(othersubset,2)-walkersubsets(othersubset,1)).*rand);   % Pick a random walker from the othersubset since 'unfortunately simultaneously advancing each walked instead of evolving walker in series subtly violates detailed balance. So we must split the full ensemble into two subsets' - Foreman-Mackey et al emcee: The MCMC Hammer
                        %otherwalker(otherwalker==ii)=otherwalker+1;                                                                             % If the randomly selected walker is the current one, choose the next walker
                        otherwalker=ceil(rand*(nwalkers-1));
                    elseif strcmp(invert.ensemble_move_style, 'walk') == 1
                        %otherwalker = ceil(walkersubsets(othersubset,1) + (walkersubsets(othersubset,2)-walkersubsets(othersubset,1)+1).*rand(3));
                        otherwalkers=ceil(rand(3,1)*(nwalkers-1));
                        %otherwalkers(otherwalkers==ii)=otherwalkers+1;
                    end

                    % Do Stretch or walk move to update current walker
                    if strcmp(invert.ensemble_move_style, 'stretch') == 1
                        K1=1/(sqrt(a)-1/sqrt(a));                                                                                                   % This is the CDF relating to their suggesting sampling equation
                        K2=1/sqrt(a);
                        z = ((rand/K1+K2).^2);                                              % Calculating a scaling which which to perform the stretch move - same value for each model parameter
                        [ m_trial(circharmparameters==0) ] = do_ensemble_stretch_move( m_curr(circharmparameters==0,:), ii, otherwalker, z); % Draw new model parameters with ensemble approach. Each row is one model parameter, each column is a walker
                    elseif strcmp(invert.ensemble_move_style, 'walk') == 1
                        [ m_trial(circharmparameters==0) ] = do_ensemble_walk_move( m_curr(circharmparameters==0,:), ii, otherwalkers);      % If solving for faultsize with circharm (model parameters 7-9) we don't want to update these every iteration
                    end
                elseif  strcmp(invert.ensemble_sampling, 'no') == 1
                    m_trial(circharmparameters==0 & m_on==1) = m_curr(circharmparameters==0 & m_on==1) + ((step_sizes(circharmparameters==0 & m_on==1) - (-step_sizes(circharmparameters==0 & m_on==1))) .*rand(sum(circharmparameters==0 & m_on==1),1) - step_sizes(circharmparameters==0 & m_on==1));
                end        % only add steps to patches that are on. m_curr remains full length, and stores old values, for when these patches turn back on
                
            end


            % Check each model parameter is within its prior     Slip=1, Alpha=2, Rake=3, Dip=4, Offset=5, Beta=6, Circharmmag=7, Circharmrotation=8, circharmxy=9, corrleationlengths=10
            toosmall = m_trial < m_min;
            m_trial(toosmall) = 2*m_min(toosmall) - m_trial(toosmall);
            toobig = m_trial > m_max;
            m_trial(toobig) =  2*m_max(toobig) - m_trial(toobig);
            
            % Sort out offset
            if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
                offset_modelparameter_temp = m_trial(m_identifyer_master==5);
                for q = 1:n_InSAR_scenes
                    offset_temp(InSAR_identifyer==q,1) = offset_modelparameter_temp(q);
                end
                offset_temp((end+1):length(d))=0;    % No offset on GPS
            else
                offset_temp = zeros( length(d), 1);
            end
                        
            % Recalculate G_ss and G_ds, if solving for dip
            if strcmp(invert.solve_for_dip, 'yes') == 1
                G_ss_temp(:, first_patch_in_strand(k):last_patch_in_strand(k)) = interp3(1:n_slip_patches_on_each_fault_strand(k), 1:n_data, dip_LUT, G_ss(:,first_patch_in_strand(k):last_patch_in_strand(k),:), 1:n_slip_patches_on_each_fault_strand(k), 1:n_data, m_trial(m_identifyer==4));
                G_ds_temp(:, first_patch_in_strand(k):last_patch_in_strand(k)) = interp3(1:n_slip_patches_on_each_fault_strand(k), 1:n_data, dip_LUT, G_ds(:,first_patch_in_strand(k):last_patch_in_strand(k),:), 1:n_slip_patches_on_each_fault_strand(k), 1:n_data, m_trial(m_identifyer==4));
            end
                        
            
            % Update slipping area, if solving for fault size
            if strcmp(invert.solve_for_fault_size, 'yes') == 1 && store_number >= faultsizechangestart && rem(store_number,faultsizechangefreq) == 0 || strcmp(invert.solve_for_fault_size, 'yes') == 1 && faultsizechangestart == 1 && store_number == 1
                
                % Calc circular harmonics with current model parameters
                [circx,circz]=circharm(m_trial(m_identifyer_master==7),m_trial(m_identifyer_master==8),0);     % Set to '1' if want to plot
                
                % Work out which patches are 'on' and 'off'
                circharm_center_trial = m_trial(m_identifyer_master==9);
                onoffidentifyer = inpolygon(patchx, patchz, (circx+circharm_center_trial(1)), (circz+circharm_center_trial(2)))';
                
                % Check with figure
                % plot_on_patches

                % Recalculate first, last, on slip patches on each fault
                % strand
                for n = 1:n_fault_strands_for_smoothing
                    n_slip_patches_ON_on_each_fault_strand_for_smoothing(n,1) = sum(onoffidentifyer(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),1));
                end
                n_slip_patches_ON_on_each_fault_strand_for_smoothing_save(i/faultsizechangefreq) = n_slip_patches_ON_on_each_fault_strand_for_smoothing;
                first_patch_in_strand_for_smoothing(1) = 0;
                first_patch_in_strand_for_smoothing(2:(n_fault_strands_for_smoothing+1),1) = cumsum(n_slip_patches_ON_on_each_fault_strand_for_smoothing);
                first_patch_in_strand_for_smoothing = first_patch_in_strand_for_smoothing+1;
                first_patch_in_strand_for_smoothing(n_fault_strands_for_smoothing+1) = [];
                last_patch_in_strand_for_smoothing = cumsum(n_slip_patches_ON_on_each_fault_strand_for_smoothing);
                
                %  Update fault size - NOTE THIS ONLY WORKS FOR ONE FAULT STRAND FOR NOW. If doing for more than one fault strands need to sum along the fault, not just the change in utmx and utmy between the patches.
                tippytop = min(disloc_model(8,onoffidentifyer==1));     % this finds the minimum top depth of any patches that are on from the matrix disloc model
                verybottom = max(disloc_model(8,onoffidentifyer==1));   % this finds the maximum top depth of any patches that are on from the matrix disloc model
                slippingpatch_width_for_smoothing =  verybottom-tippytop;
                
                farleft = min(disloc_model(1,onoffidentifyer==1));
                farright = max(disloc_model(1,onoffidentifyer==1));
                slippingpatch_length_for_smoothing = abs(farleft-farright);           % abs in case in a synthetic test you've defined +x and -x
                
                if length(slippingpatch_width_for_smoothing) == 0
                    slippingpatch_width_for_smoothing = 0;
                end
                
                if length(slippingpatch_length_for_smoothing) == 0
                    slippingpatch_length_for_smoothing = 0;
                end
                
                if strcmp(predominant_faulting_style, 'ss') == 1
                    % Strike slip parameters
                    a_as_ss_mean = 1860 + 0.34*slippingpatch_length_for_smoothing;
                    a_as_ss_std = sqrt(1120^2+0.03^2*slippingpatch_length_for_smoothing^2);          % Errors from (Mai and Beroza, 2002), Table 2, Pg 14
                    a_dd_ss_mean = -390 + 0.44 *slippingpatch_width_for_smoothing;                   % Propagation of errors from talbe on Wikipedia (shh)  https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                    a_dd_ss_std = sqrt(470^2+0.04^2*slippingpatch_width_for_smoothing^2);
                elseif strcmp(predominant_faulting_style, 'ds') == 1
                    % Dip slip parameters
                    a_as_ds_mean =  1100 + 0.31*slippingpatch_length_for_smoothing;
                    a_as_ds_std = sqrt(400^2+0.01^2*slippingpatch_length_for_smoothing^2);
                    a_dd_ds_mean =  580 + 0.35*slippingpatch_width_for_smoothing;
                    a_dd_ds_std = sqrt(240^2+0.01^2*slippingpatch_width_for_smoothing^2);
                end
                
                % Recalculate correlation lengths
                for n = 1 : n_fault_strands_for_smoothing
                    if strcmp(predominant_faulting_style(n), 'ss') == 1
                        if strcmp(invert.solve_for_correlation_length, 'yes') == 1
                            a_as(n) = a_as_ss_mean + erfinv(m_trial(m_identifyer_master==10))*sqrt(2)*a_as_ss_std;   % Have to use error function trick to draw from a normal distribution, since we can't do a random walk with randn function
                            a_dd(n) = a_dd_ss_mean + erfinv(m_trial(m_identifyer_master==11))*sqrt(2)*a_dd_ss_std;
                        else
                            a_as(n) =  1860 + 0.34*(slippingpatch_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
                            a_dd(n) =  -390 + 0.44 * (slippingpatch_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002)
                        end
                    elseif strcmp(predominant_faulting_style(n), 'ds') == 1
                        if strcmp(invert.solve_for_correlation_length, 'yes') == 1
                            a_as(n) = a_as_ds_mean + erfinv(m_trial(m_identifyer_master==10))*sqrt(2)*a_as_ds_std;
                            a_dd(n) = a_dd_ds_mean + erfinv(m_trial(m_identifyer_master==11))*sqrt(2)*a_dd_ds_std;
                        else
                            a_as(n) =  1100 + 0.31*(slippingpatch_length_for_smoothing(n));          % Along-strike.  fault_length = meters  (Mai and Beroza, 2002)
                            a_dd(n) =  580 + 0.35* (slippingpatch_width_for_smoothing(n));         % Down-dip.      fault_width  = meters  (Mai and Beroza, 2002)
                        end
                    end
                end

                % Recalculate scaled distances
                [r_over_a, ~] = calc_scaled_dist( n_fault_strands_for_smoothing, disloc_model, a_as, a_dd, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master, along_strike_sep_dist, n_along_strike_patches_for_smoothing, n_down_dip_patches_for_smoothing,fault_strand_togetherness);

                % Recalculate sigma_s
                sigma_s = calc_sigma_s( r_over_a(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)), H(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)));       % NOTE that this is actually a SEPARATE SIGMA S MATRIX FOR EACH FAULT STRAND, just stored in one big matrix
%                %Uncomment below for multiple fault strands 
%                 sigma_s_master = [];        % sigma_s_master is updated because faultlength and faultwidth changed, so we had to recalculate everything
%                 for n = 1: n_fault_strands_for_smoothing
%                     sigma_s_temp = calc_sigma_s( r_over_a(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)), H(first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n),first_patch_in_strand_for_smoothing_master(n):last_patch_in_strand_for_smoothing_master(n)));       % NOTE that this is actually a SEPARATE SIGMA S MATRIX FOR EACH FAULT STRAND, just stored in one big matrix
%                     sigma_s_master = blkdiag(sigma_s_master,sigma_s_temp);
%                     clear sigma_s_temp
%                 end
%                 sigma_s = sigma_s_master;   % sigma_s only contains the rows and columns that correspond to patches that are 'on'.   
                
                % add correlation matrix stabiliser, if necessary
                if strcmp(invert.add_correlation_matrix_stabiliser, 'yes') == 1
                    sigma_s = sigma_s + diag( 0.01*ones(total_n_slip_patches,1));
                end
                sigma_s(:, onoffidentifyer==0) = [];            % slice columns out of master
                sigma_s(onoffidentifyer==0, :) = [];            % slice rows out of master

                det_sigma_s = det(sigma_s);
                inv_sigma_s = inv(sigma_s);     %Uncomment below for multiple fault strands   
                %det_sigma_s = [];
                %inv_sigma_s = [];
%                 for n = 1: n_fault_strands_for_smoothing            % NOTE that this is actually a SEPARATE invSIGMA S MATRIX FOR EACH FAULT PATCH, just stored in one big matrix, in the third dimension for each fault strand
%                     det_sigma_s(n,1) = det(sigma_s(first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n),first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n)));
%                     inv_sigma_s_temp(:, :) = inv(sigma_s(first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n),first_patch_in_strand_for_smoothing(n):last_patch_in_strand_for_smoothing(n)));
%                     inv_sigma_s = blkdiag(inv_sigma_s, inv_sigma_s_temp);
%                     inv_sigma_s_temp = [];
%                 end
                clear inv_sigma_s_temp
                
                det_sigma_s_keep(:,faultsizecount) = det_sigma_s;
                det_sigma_s_ratio = prod((det_sigma_s_keep(:,faultsizecount)/det_sigma_s_keep(:,faultsizecount-1))^(-0.5));   % But we only need to calculate this for the change between n patches being on and m patches being on. once a trial is accepted, m patches will remain on for the next 1000 iterations. so once a trial is accepted we need to set this back to 1
                m_on = [onoffidentifyer; ones(n_fault_strands_for_smoothing,1); onoffidentifyer; ones(size(dip_initial,1),1); ones(size(offset_modelparameter_initial,1),1); ones(size(beta_initial,1),1); ones(sum(circharmparameters), 1)];
                m_identifyer = m_identifyer_master;
                m_identifyer(m_on==0) = [];  % Remove appropriate slip and rake rows from m_identifyer
                m_on_keep(:,faultsizecount) = m_on;
                
                
%               % Look at features of newly turned on/off patches
                births = zeros(total_n_slip_patches,1);
                births(m_on_keep(m_identifyer==1,faultsizecount)==1 & m_on_keep(m_identifyer==1,faultsizecount-1)==0) = 1;      % slip patches that are on for the current faultsizecount, but were off for the one before
                deaths = zeros(total_n_slip_patches,1);
                deaths(m_on_keep(m_identifyer==1,faultsizecount)==0 & m_on_keep(m_identifyer==1,faultsizecount-1)==1) = 1;    % slip patches that are off for the current faultsizecount, but were on for the one before

                % Calculate proposal ratio
                %proposal_ratio = ((total_n_slip_patches-n_slip_patches_ON_on_each_fault_strand_for_smoothing)/(n_slip_patches_ON_on_each_fault_strand_for_smoothing+1))^sum(births) * ((n_slip_patches_ON_on_each_fault_strand_for_smoothing)/(total_n_slip_patches+ n_slip_patches_ON_on_each_fault_strand_for_smoothing+1))^sum(deaths)

                % Calculte the pre_prior_ratio
                transdimensional_prior_ratio = 1 / ((slip_range)^(sum(births)-sum(deaths))*(rake_range)^(sum(births)-sum(deaths)));
            end
                
                % Recalculate G matrix, if solving for rake or dip
                if strcmp(invert.variable_rake, 'yes')==1 || strcmp(invert.solve_for_dip, 'yes') == 1
                    G_temp = G_ss_temp(:, onoffidentifyer==1)*diag(cosd(m_trial(m_identifyer_master==3&m_on==1))) +  G_ds_temp(:, onoffidentifyer==1)*diag(sind(m_trial(m_identifyer_master==3&m_on==1)));      % just for on rake patches
%                     if i == sens_test(sens_test_number) && ii == 1     % If this is a sensitivity iteration then save this G matrix for sensitivity calculation later
%                         G_sensitivity_test = G_temp;
%                     end
                end
                
                % Remove 'off' patches, and corresponding rakes, and G_temp
                if sum(onoffidentifyer) ~= total_n_slip_patches
                    m_trial(m_on==0) = [];
                end
                
                
                
                % Calculate moment and moment likelihood
                if strcmp(invert.regularise_moment, 'yes') == 1
                    M0_temp = sum((elastic_params.mu_okada) * spatial_model2column(onoffidentifyer==1) .* spatial_model3column(onoffidentifyer==1) .* (m_trial(m_identifyer==1)));
                    M0_likelihood_temp = normpdf( M0_temp, data.seismic_moment, data.moment_std);
                end
                
                % Beta ratio
                if strcmp(invert.solve_for_beta, 'yes') == 1
                    beta_ratio = (m_trial(m_identifyer==6)^2 / m_curr(m_identifyer_master==6)^2) ^ (-length(d)/2);     % Else this is set at 1 elsewhere
                end
                
                % SMOOTHED
                if strcmp(invert.inversion_type, 'bayesian') == 1 && strcmp(invert.smoothing, 'VK') || strcmp(invert.smoothing, 'laplacian') || strcmp(invert.smoothing, 'minimumnorm') == 1
                    
                    % Sort out alpha^2
                    %for l = 1: n_fault_strands_for_smoothing
                    if strcmp(priors.alpha2_prior, 'logarithmic') == 1
                        alpha2_trial = 10.^m_trial(m_identifyer==2);             % m_identifyer has 'off' values removed
                    else
                        alpha2_trial = m_trial(m_identifyer==2);
                    end
                    %end
                    alpha2_ratio = prod( ((alpha2_trial.^(n_slip_patches_ON_on_each_fault_strand_for_smoothing) ./ alpha2_curr.^(previous_n_slip_patches_ON_on_each_fault_strand_for_smoothing)) .^ (-1/2)));      % do the product because we need to multiply all the ratios for the n fault strands
                    
                    % Calculate likelihood
                    if strcmp(invert.smoothing, 'VK') == 1    % VK ***************
                        [logprior_trial, singularflag] = calc_logprior_VK( m_trial(m_identifyer==1), inv_sigma_s, alpha2_trial, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing); % Note: this is just calculating the exponent
                        exponent_ratio = prod ( exp( logprior_trial - logprior_curr(:,ii)) );                                    % do the product because we need to multiply all the ratios for the n fault strands
                        prior_ratio = ( alpha2_ratio * exponent_ratio *det_sigma_s_ratio );      % det_sigma_s_ratio is set at 1 before the iterations between, and is only updated if solving for fault size (and therefore det_sigma_s changes between iterations)
                    elseif strcmp(invert.smoothing, 'laplacian') == 1 % Laplacian****************
                        logprior_trial = calc_logprior_laplacian( m_trial(m_identifyer==1), L, alpha2_trial, n_fault_strands_for_smoothing, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master);
                        exponent_ratio = prod ( exp( logprior_trial - logprior_curr(:,ii)) );         % do the product because we need to multiply all the ratios for the n fault strands
                        prior_ratio = ( alpha2_ratio * exponent_ratio );
                    elseif strcmp(invert.smoothing, 'minimumnorm') == 1 % Minmum norm ***************
                        for q = 1:n_fault_strands_for_smoothing
                            s_temp = slip_temp(first_patch_in_strand_for_smoothing(q):last_patch_in_strand_for_smoothing(q));
                            logprior_trial(q,1) =  ( -1/(2*alpha2_trial(q)) * sum(s_temp'*s_temp));
                        end
                        exponent_ratio = prod(exp( logprior_trial - logprior_curr(:,ii)));
                        prior_ratio = ( alpha2_ratio * exponent_ratio );
                    end
                    
                    %                 if singularflag==1
                    %                     keyboard
                    %                 end
                    
                    prior_ratio = prior_ratio * transdimensional_prior_ratio;
                    
                    % Calculate prior ratio
                    if strcmp(invert.regularise_moment, 'yes') ==1
                        M0_likelihood_ratio =  M0_likelihood_temp / M0_likelihood_curr;
                        prior_ratio = prior_ratio * M0_likelihood_ratio;
                    end
                    
                    if strcmp(invert.ensemble_move_style, 'stretch') == 1
                        prior_ratio = z^(n_model_parameters-1)*prior_ratio;
                    end
                    
                    %                 if strcmp(invert.solve_for_fault_size, 'yes') == 1
                    %                     prior_ratio = likelihood_ratio * proposal_ratio;
                    %                 end
                    
                % UNSMOOTHED
                elseif strcmp(invert.smoothing, 'none') == 1
                    prior_ratio = 1;    % Always want to accept the prior, if no priors
                    if strcmp(invert.regularise_moment, 'yes') ==1
                        M0_likelihood_ratio = M0_likelihood_temp / M0_likelihood_curr;
                        prior_ratio = prior_ratio * M0_likelihood_ratio;
                    end
                end
                
                
                % MCMC Step 1 - Sampling the prior
                if prior_ratio > rand    % Take metropolis step. Or not, as the case may be
                    
                    true_trials = true_trials+1;
                    %prior_ratio_keep(store_number+1, 1) = prior_ratio;
                    
                    % Beta- if not solving for it beta_ratio stays at 1
                    if strcmp(invert.solve_for_beta, 'yes') ==1
                        beta_ratio = (m_trial(m_identifyer==6)^2 / m_curr(m_identifyer_master==6)^2) ^ (-length(d)/2);
                    end
                    
                    % Calculate the likelihood
                    logL_trial = calc_loglikely(m_trial(m_identifyer==1), d, G_temp, inv_sigma_d, offset_temp, m_trial(m_identifyer==6) );
                    likelihood_ratio = beta_ratio * exp( logL_trial - logL_curr(:,ii));
                    
                    % MCMC Step 2 - Sampling the posterior (using the likelihood)
                    if  likelihood_ratio > rand
                        
                        % Update
                        m_curr(m_on==1,ii) = m_trial;
                        logprior_curr(:,ii) = logprior_trial;
                        logL_curr(:,ii) = logL_trial;
                        det_sigma_s_ratio = 1;  % Whilst we were changing the number of 'on' fault patches det_sigma_s is different. But once we've changed fault patches and accepted a trial then the number of fault patches is the same between iterations. So this is back to 1. And if we're not changing the number of fault patches then this is always right
                        alpha2_curr = alpha2_trial;
                        if strcmp(invert.regularise_moment, 'yes') ==1
                            M0_curr = M0_temp;
                            M0_likelihood_curr = M0_likelihood_temp;
                        end
                        if strcmp(invert.solve_for_dip, 'yes') ==1
                            G_ss_curr = G_ss_temp;
                            G_ds_curr = G_ds_temp;
                        end
                        if strcmp(invert.smoothing, 'VK') == 1
                            logposterior_curr =  sum( log((2*pi*alpha2_curr).^(-n_slip_patches_ON_on_each_fault_strand_for_smoothing/2))) + sum(logprior_curr(:,ii)) + logL_curr;
                        elseif strcmp(invert.smoothing, 'none') == 1
                            logposterior_curr = logL_curr(:,ii);
                        elseif strcmp(invert.smoothing, 'VK') ~= 1
                            logposterior_curr =  sum(logprior_curr(:,ii)) + logL_curr;
                        end
                        if strcmp(invert.regularise_moment, 'yes') ==1
                            logposterior_curr =  logposterior_curr + log(M0_likelihood_curr);
                        end
                        if strcmp(invert.solve_for_fault_size, 'yes') == 1
                            mcircharmfreeze = m_curr(circharmparameters==1);
                            previous_n_slip_patches_ON_on_each_fault_strand_for_smoothing = n_slip_patches_ON_on_each_fault_strand_for_smoothing;
                            %transdimensional_prior_ratio = 1;  % once a new transdimensional trial has been accepted, the pre_prior_ratio is 1 till the next change.
                            %proposal_ratio = 1;  % once a new transdimensional trial has been accepted, the proposal ratio is 1 till the next change.
                        end
                        
                        %if strcmp(invert.slip_prior, 'gaussian') + strcmp(invert.slip_prior, 'logarithmic')  > 0
                        %            b_curr = b_trial;                                  % Update param_curr, since this is the parameter that we're actually drawing, so we need to do our random walk from this
                        % end
                        
                    else % Do not use this temporary slip distribution as a trial; calculate new trial.
                        likelihood_rejections(1,store_number) = 1;                             % trial_rejections = 1 if trial rejected, remains 0 if trial was accepted (but this does not necessarily mean that the step will be accepted...).
                    end
                    
                else
                    prior_rejections(1,store_number) =  1;              % trial_rejections = 1 if step rejected, remains 0 otherwise. NOT 0 DOES NOT NECESSARILY EQUAL ACCEPTANCE, since it might remain 0 because the trial was rejected, before we even got to this step. it wasn't given the chance to be accepted or rejected. i can explain this better but i have to go to badminton.
                end
                
                % Store. Either the updated trial or the previous trial. We store the 'curr' values - so this will either be the updated curr if a trial was accepted, or it will be the same as the previous iteration if the trial was rejected.
                if ii == nwalkers
                   m_keep(m_on==1,i,:) = m_curr(m_on==1);      % On patches saved, the rest becomes 0 (both slip and rake)
                   %m_keep(m_on==0,i,:)=0;
                end
                if strcmp(invert.solve_for_fault_size, 'yes') == 1
                    m_trial = m_curr(:,ii); % This is to set m_trial back to the size of total_n_model_parameters if it's changed whilst solving for fault size
                end
                if i > burn_in_remove_number
                    onoffidentifyerkeep = onoffidentifyerkeep + onoffidentifyer;
                end
                if strcmp(invert.regularise_moment, 'yes') == 1
                    M0_keep(1,store_number) = M0_curr;
                    M0_likelihood_keep(1,store_number) = M0_likelihood_curr;
                end
                logprior_keep(:,store_number) = sum(log((2*pi*alpha2_curr).^(-n_slip_patches_ON_on_each_fault_strand_for_smoothing/2))) + sum(logprior_curr(:,ii));   % NOTE: not including det_sigma_s
                logL_keep(1,store_number) = logL_curr(:,ii);                                                            % NOTE: not including beta
                logposterior_keep(1,store_number) = logposterior_curr(:,ii);  % One column for each iteration, containing posterior for each patch (number of rows = number of slip patches)

                transdimensional_prior_ratio = 1;                % At the end of a transdimensional perturbation iteration, update the ratio so that the model parameter perturbations aren't rejected due to this ratio.
     


            % .........................................................................  

            % *********************************************************************** %
            % ......................................................................... 


            % Announce progress

                if i == abs(iterations/4) && ii == 1    % So that it only tells us this for the first walker
                    disp('    ')
                    disp('One quarter of the way through...')
                    disp('..')
                    disp('    ')
                elseif i == abs(iterations/2) && ii == 1
                    disp('    ')
                    disp('Halfway through...')
                    disp('...')        
                    disp('    ') 
                elseif i == abs(iterations*3/4) && ii == 1
                    disp('    ')
                    disp('Three quarters through...')     
                    disp('....')
                    disp('    ')
                end


            % ............................................................. 

            % % Sensitivity test ******************************************
            % .............................................................     

            if i == sens_test(sens_test_number) && ii == 1 % So we only do this for the first walker in iteration 100. Not for all walkers in iteration 100.
                
                % Calculate total rejections
                total_rejections = likelihood_rejections + prior_rejections + pre_prior_rejections;      % This is if the random step was rejected at any point - prior ratio or likelihood ratio.
                store_number_at_sens_test(1,sens_test_number) = store_number;
                
                % Calculate how many patches are currently 'on'
                if strcmp(invert.solve_for_fault_size, 'yes') == 1
                    disp(' ')
                    disp(['Total number of slip patches on = ', num2str(sum(onoffidentifyer))])
                end
                
                % Work out rejection rate - for ensemble sampling
                if strcmp(invert.ensemble_sampling, 'yes') == 1
                    if sens_test_number == 1
                        disp(['Current rejection rate is ', num2str( sum(total_rejections(1:store_number)) / store_number), ' with ', num2str(store_number), ' trials'  ])
                    else
                        disp(['Current rejection rate is ', num2str( sum(total_rejections(store_number_at_sens_test(sens_test_number-1):store_number)) / (store_number_at_sens_test(sens_test_number)-store_number_at_sens_test(sens_test_number-1))), ' with ', num2str(store_number), ' trials'  ])
                    end
                    
                else
                    
                    % Work out rejection rate - not for ensemble sampling
                    i_at_sens_test(1,sens_test_number) = i;                % because an interation is only at true iteration after each walker has been perturbed. so i will equal 100 for nwalker times before one real iteration has been completed
                    
                    if sens_test_number == 1            % on the first sensitivity test only
                        rej_ratio_since_last_sens_test =  (sum(total_rejections(1:(sens_test(sens_test_number)-1))))  / (sens_test(sens_test_number) );  % Check the rejection ratio for the iterations since the last sensitivity test. This is the number of rejections in the last batch of iterations, divided by the number of iterations.
                    else
                        rej_ratio_since_last_sens_test =  sum(total_rejections(sens_test(sens_test_number-1):(sens_test(sens_test_number))-1)) / (sens_test(sens_test_number) - sens_test(sens_test_number-1) );  % the first minus one is because
                    end
                    
                    if strcmp(invert.smoothing, 'none') == 0
                        posterior_accepts_at_sens_test = true_trials - sum(likelihood_rejections);   % Only a true trial if passed prior loop
                    else
                        accepts_at_sens_test = i - sum(total_rejections);
                    end
                    
                    if strcmp(invert.smoothing, 'none') == 0
                        disp(['Current rejection rate = ' num2str(rej_ratio_since_last_sens_test), ' with ' num2str(posterior_accepts_at_sens_test) ' posterior accepts after ', num2str(true_trials), ' true trials (' num2str(i), ' iterations)']);
                    elseif  strcmp(invert.smoothing, 'none') == 1
                        disp(['Current rejection rate = ' num2str(rej_ratio_since_last_sens_test), ' with ' num2str(accepts_at_sens_test) ' accepts after ', num2str(i), ' iterations']);
                    end
                    
                    % Update probability target using current rejection ratio
                    rej_ratio_save(:,sens_test_number) = rej_ratio_since_last_sens_test;
                    probability_target = probability_target / (rej_ratio_since_last_sens_test / rejection_target); % Compare rejection ratio to ideal rejection ratio and make ideal_probability perturbation bigger/smaller, depending on if rej_ratio is too low/high
                    probability_target(probability_target>1) = 0.99;        % If this gets to be > 1 this messes up step size calculation like crazy.
                    
                    G_sensitivity_test = G_ss_temp(:, onoffidentifyer==1)*diag(cosd(m_trial(m_identifyer_master==3&m_on==1))) +  G_ds_temp(:, onoffidentifyer==1)*diag(sind(m_trial(m_identifyer_master==3&m_on==1)));      % just for on rake patches
                    
                    % Calculate new step sizes
                    if singularflag == 0 & i < 10000           % only do step sizes test if current matrix isn't singular, if changing faultsize
                        [new_step_sizes] = calculate_step_sizes( m_curr, total_n_slip_patches, invert, G_sensitivity_test, d, inv_sigma_d, inv_sigma_s_master, step_sizes, probability_target, G_ss_curr, G_ds_curr, n_fault_strands, n_fault_strands_for_smoothing, n_slip_patches_on_each_fault_strand, first_patch_in_strand, last_patch_in_strand, first_patch_in_strand_for_smoothing_master, last_patch_in_strand_for_smoothing_master, det_sigma_s, L, M0_likelihood_curr, elastic_params, spatial_model2, spatial_model3, data , m_identifyer_master, G_ss, G_ds, dip_LUT, n_data, fault_strand_identifyer, priors.min_dip, priors.max_dip, onoffidentifyer, n_slip_patches_ON_on_each_fault_strand_for_smoothing, n_slip_patches_on_each_fault_strand_for_smoothing, d_identifyer, patchx, patchz, disloc_model, n_down_dip_patches_for_smoothing, n_along_strike_patches_for_smoothing, along_strike_sep_dist, n_along_strike_patches, n_down_dip_patches, fault_strand_togetherness, H,first_patch_in_strand_for_smoothing,last_patch_in_strand_for_smoothing, fault_length_for_smoothing, fault_width_for_smoothing, InSAR_identifyer, n_harmonics, predominant_faulting_style, priors, n_InSAR_scenes, m_on, circharmparameters);
                    end
                    step_sizes_keep(:, sens_test_number) = new_step_sizes;
                    step_sizes = new_step_sizes;
                    
                end
                
                sens_test_number = sens_test_number + 1;    % Ready for the next sensitivity test.
                
            end
        end
        
        
        
    end



    %% *************************************************************************
    % Enough looping now ****************************************************** 
    % ************************************************************************* 

    % Do some bayesian tidying up .............................................
    disp('  ');
    disp('Almost done. Tidy tidy tidy.');


    % Remove the burn in. I don't have a sophisticated way to find the number of steps to remove. So I just remove the first 2000
    if strcmp(invert.load_old_MCMC_chain, 'yes') == 0  % Don't need to remove burn-in if you're loading an old chain, presumably.
        
        %burn_in_remove_number = 1:(burn_in_remove_number);
        m_keep(:,1:burn_in_remove_number/nwalkers,:)=[];
        logposterior_keep(:,1:burn_in_remove_number) = [];
        if strcmp(invert.regularise_moment, 'yes') == 1
                M0_keep(:,1:burn_in_remove_number) = [];
                M0_likelihood_keep(:,1:burn_in_remove_number) = [];
        end
        logL_keep(:,1:burn_in_remove_number)=[];
        
    end  
    
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
        
        % Expand matrices if doing ensemble_sampling
        if strcmp(invert.ensemble_sampling, 'yes') == 1
            slip_keep = reshape(slip_keep, [total_n_slip_patches, (size(slip_keep,2)*size(slip_keep,3))]);
            alpha2_keep = reshape(alpha2_keep, [n_fault_strands_for_smoothing, (size(alpha2_keep,2)*size(alpha2_keep,3))]);
            rake_keep = reshape(rake_keep, [total_n_slip_patches, (size(rake_keep,2)*size(rake_keep,3))]);
            if strcmp(invert.solve_for_dip, 'yes') ==1
                dip_keep = reshape(dip_keep, [total_n_slip_patches, (size(dip_keep,2)*size(dip_keep,3))]);
            end
            if strcmp(invert.solve_for_beta, 'yes') == 1
                beta_keep = reshape(beta_keep, [1, (size(beta_keep,2)*size(beta_keep,3))]);
            end
            if strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
                offset_keep = reshape(offset_keep, [n_InSAR_scenes, (size(offset_keep,2)*size(offset_keep,3))]);
            end
            if strcmp(invert.solve_for_fault_size, 'yes') == 1;
                circharm_coeffs_keep = reshape(circharm_coeffs_keep, [nharmonics, (size(circharm_coeffs_keep,2)*size(circharm_coeffs_keep,3))]);
                circharm_phi_keep = reshape(circharm_phi_keep, [nharmonics-1, (size(circharm_phi_keep,2)*size(circharm_phi_keep,3))]);
                circharm_center_keep = reshape(circharm_center_keep, [2, (size(circharm_center_keep,2)*size(circharm_center_keep,3))]);
            end
        end
end    % End of bayesian inversion


%% After the inversions

disp('  ');
disp('.....Okay, enough of that now.');
disp('  ');

inversion_time = toc;

%% Save everything - in case display_result doesn't run proper

if strcmp(invert.pad_edges_with_zeros, 'yes') == 1
     savename = [housekeeping.save_name, '_', num2str(n_down_dip_patches_for_smoothing(1)), 'x', num2str(n_along_strike_patches_for_smoothing(1)), '_', invert.smoothing, 'smooth_', invert.solve_for_dip, 'dip_', num2str(invert.iterations), '_', invert.regularise_moment, 'M0reg_', priors.slip_prior, '_zeropadded.mat'];
else
     savename = [housekeeping.save_name, '_', num2str(n_down_dip_patches_for_smoothing(1)), 'x', num2str(n_along_strike_patches_for_smoothing(1)), '_', invert.smoothing, 'smooth_', invert.solve_for_dip, 'dip_', num2str(invert.iterations), '_', invert.regularise_moment, 'M0reg_', priors.slip_prior, '_', invert.solve_for_fault_size, 'patchesonoff'];            
end

%save(savename);
save(savename, '-v7.3');


%% Calculate your averages, display the result

if strcmp(invert.inversion_type, 'least_squares') == 1 && strcmp(invert.smoothing, 'tikhonov') == 1
%         visual_least_sq_solution = reshape(least_sq_solution, n_down_dip_patches, n_along_strike_patches);
%         imagesc(visual_least_sq_solution);
%         xlabel('Along Strike (km)'); ylabel('Down dip (km)'); ylabel(colorbar, 'Slip (m)');        
%             if strcmp(invert.smoothing, 'none') == 1;
%                             smoothstring = ('no');
%             elseif strcmp(invert.smoothing, 'laplacian') == 1;
%                 smoothstring = ('Laplacian');
%             elseif strcmp(invert.smoothing, 'VK') == 1;
%                 smoothstring = ('von Karman');
%             end    
%         title(['Least-squares solution (with ', smoothstring, ' smoothing)'])
%         axis equal tight
%         colorbar
        display.plotsurfacedisp = 'yes';
        G_mostlikely = best_G_tik;
        patch_mode = best_slip_tik;
        patch_MAP = best_slip_tik;
        
        fontsize_plot = 22;
        set(0,'DefaultAxesFontSize',fontsize_plot)
        set(0,'defaultAxesFontName', 'Times New Roman')
        display.calc_confidence = 'no';
        display.plotmean = 'no';
        display.plotmode = 'no';
        display.plotmedian = 'no';
        display.plotmaxlikely = 'no';
        display.plotallslips = 'no';
        display.plotprob = 'no';
        display.plothists = 'no';
        display.plotmarginalPDFs = 'no';
        display.plotMAP = 'no';
        display.plot_resolution_matrix = 'no';
        display.plot3d = 'no';
        run('display_result.m');
        
        
        figure('position', [100, 350, 1600, 1200])
        faults = disloc_model; 
        faults(6,:) = best_slip_tik'; 
        doplot3d(faults', 'jet');
        hold on;
        colorbar
        title('3D mode')
        xlabel('UTM x')
        ylabel('UTM y')
        zlabel('Depth (km)')
        title('tikhonov')
        % Add rake
        quiv_mags = [patch_mode.*cosd(best_rake_tik), patch_mode.*sind(best_rake_tik)];
        [xquivmag,yquivmag,zquivmag] = reproject_quiv(quiv_mags,disloc_model(3,:)',disloc_model(4,:)');   % nicked from Tom.   reproject_quiv(ss_mag, ds_mag, strike, dip)
        scale_factor = 1;
        quiver3(disloc_model(1,:)'/1000,disloc_model(2,:)'/1000,-0.5*(disloc_model(8,:)+disloc_model(9,:))'/1000,xquivmag*scale_factor, yquivmag*scale_factor, zquivmag*scale_factor, 'k', 'Linewidth', 1.5, 'Autoscale', 'off');
    
        
        cutoff=0.95*sum(best_slip_tik);
        [I]= find(cumsum(sort(best_slip_tik, 'descend'))>cutoff,1, 'first');
        slipcutoff=sort(best_slip_tik(I));
        
elseif strcmp(invert.inversion_type, 'SVD') == 1
    
        figure('position', [100, 350, 1600, 1200])
        faults = disloc_model; 
        faults(6,:) = slip_svd'; 
        doplot3d(faults', 'jet');
        hold on;
        colorbar
        title('3D mode')
        xlabel('UTM x')
        ylabel('UTM y')
        zlabel('Depth (km)')
        title('SVD')
        % Add rake
        quiv_mags = [slip_svd.*cosd(rake_svd), slip_svd.*sind(rake_svd)];
        [xquivmag,yquivmag,zquivmag] = reproject_quiv(quiv_mags,disloc_model(3,:)',disloc_model(4,:)');   % nicked from Tom.   reproject_quiv(ss_mag, ds_mag, strike, dip)
        scale_factor = 0.7;
        quiver3(disloc_model(1,:)'/1000,disloc_model(2,:)'/1000,-0.5*(disloc_model(8,:)+disloc_model(9,:))'/1000,xquivmag*scale_factor, yquivmag*scale_factor, zquivmag*scale_factor, 'k', 'Linewidth', 1.5, 'Autoscale', 'off');
    
        display.plotsurfacedisp = 'yes';
        G_mostlikely = G_svd;
        patch_mode = [ss_slip_svd; ds_slip_svd];
        
elseif strcmp(invert.inversion_type, 'least_squares') == 1 && strcmp(invert.smoothing, 'cv') == 1
        
        display.plotsurfacedisp = 'yes';
        G_mostlikely = G;
        patch_mode = best_slip_cv;
        patch_MAP = best_slip_cv;
        
        fontsize_plot = 22;
        set(0,'DefaultAxesFontSize',fontsize_plot)
        set(0,'defaultAxesFontName', 'Times New Roman')
        display.calc_confidence = 'no';
        display.plotmean = 'no';
        display.plotmode = 'no';
        display.plotmedian = 'no';
        display.plotmaxlikely = 'no';
        display.plotallslips = 'no';
        display.plotprob = 'no';
        display.plothists = 'no';
        display.plotmarginalPDFs = 'no';
        display.plotMAP = 'no';
        display.plot_resolution_matrix = 'no';
        display.plot3d = 'no';
        run('display_result.m');
        
        
        figure('position', [100, 350, 1600, 1200])
        faults = disloc_model; 
        faults(6,:) = best_slip_cv'; 
        doplot3d(faults', 'jet');
        hold on;
        colorbar
        title('3D mode')
        xlabel('UTM x')
        ylabel('UTM y')
        zlabel('Depth (km)')
        title('least-squares solution')
        % Add rake
        quiv_mags = [patch_mode.*cosd(rake_mode), patch_mode.*sind(rake_mode)];
        [xquivmag,yquivmag,zquivmag] = reproject_quiv(quiv_mags,disloc_model(3,:)',disloc_model(4,:)');   % nicked from Tom.   reproject_quiv(ss_mag, ds_mag, strike, dip)
        %[xquivmag,yquivmag,zquivmag] = pol2cart(disloc_model(3,:)',rake_mode,patch_mode);
        scale_factor = 0.7;
        quiver3(disloc_model(1,:)'/1000,disloc_model(2,:)'/1000,-0.5*(disloc_model(8,:)+disloc_model(9,:))'/1000,xquivmag*scale_factor, yquivmag*scale_factor, zquivmag*scale_factor, 'k', 'Linewidth', 1.5, 'Autoscale', 'off');

else
        disp('Displaying result.');
        run('display_result.m');
end

disp(['slipBERI completed in ', num2str(toc),' seconds. You''re welcome.']);
disp('  ');

p = profile('info'); % if you wanna view the profile, do >> profview(0,p)

%% Save everything again now you've calculated more things

if strcmp(invert.pad_edges_with_zeros, 'yes') == 1
     savename = [housekeeping.save_name, '_', num2str(n_down_dip_patches_for_smoothing(1)), 'x', num2str(n_along_strike_patches_for_smoothing(1)), '_', invert.smoothing, 'smooth_', invert.solve_for_dip, 'dip_', num2str(invert.iterations), '_', invert.regularise_moment, 'M0reg_', priors.slip_prior, '_zeropadded.mat'];
else
     savename = [housekeeping.save_name, '_', num2str(n_down_dip_patches_for_smoothing(1)), 'x', num2str(n_along_strike_patches_for_smoothing(1)), '_', invert.smoothing, 'smooth_', invert.solve_for_dip, 'dip_', num2str(invert.iterations), '_', invert.regularise_moment, 'M0reg_', priors.slip_prior, '_', invert.solve_for_fault_size, 'patchesonoff'];            
end

%save(savename);
save(savename, '-v7.3');

%% Now, the end is near... ************************************************

disp('  ');
disp('And that''s that.');

% A parting message - idea from Funning 2005
a=round(rand*10);
switch a
	case{0}
		disp('Parting is such sweet sorrow...')
        disp('that I’ll say good night until tonight becomes tomorrow.')
	case{1}
		disp('I''d kill for a custard cream.')
	case{2}
		disp('Not sure I''ve slipped enough puns in.');
	case{3}
		disp('Let''s do that AGAIN.');
	case{4}
		disp('By day, inverting quakes. By night, baking cakes.')
	case{5}
		disp('This whole business is the start of a slippery slope...');
	case{6}
		disp('GOOD JOB.');
	case{7}
		disp('For what do we live, but to make sport for our neighbours,');
        disp('and laugh at them in our turn?');
	case{8}
		disp('Did you hear? A lorry-load of tortoises crashed into a trainload of terrapins. What a turtle disaster.');
	case{9}
		disp('I''d rather be walking the dales.')
	case{10}
		disp('And with that, I''ll be off then.');
end

disp('   ')
disp('Keyboard mode now. To terminate keyboard mode and end the slipBERI function, type ''dbcont'' and press Enter')
keyboard

 
end

