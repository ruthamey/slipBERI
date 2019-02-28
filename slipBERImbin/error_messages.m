% Script to check all inputs for slipBERI are fine and dandy - I'll shout
% at you if they're not
%
%  18-aug-2014  rmja  Moved from slipBERI body to separate script, tidytidy

disp('  ')  
disp('Just doing some quick checks...')

if isstruct(data) + isstruct(display) +  isstruct(elastic_params) + isstruct(fault) + isstruct(housekeeping) + isstruct(invert) + isstruct(testing) ~= 7    % if one of these is true, isstruct = 1
    error('ERROR - one of your input structures isn''t a structure. I can''t help if you if you don''t help me.');
end

if strcmp(invert.inversion_type, 'least_squares') + strcmp(invert.inversion_type, 'bayesian')+ strcmp(invert.inversion_type, 'SVD') < 1        % if one of these is true, strcmp = 1
    error('ERROR - inversion_type must either be ''least_squares'' or ''bayesian'' or ''SVD''. Let''s try again, shall we?');
end


if strcmp(invert.smoothing, 'none') + strcmp(invert.smoothing, 'laplacian') + strcmp(invert.smoothing, 'VK') + strcmp(invert.smoothing, 'minimumnorm') + strcmp(invert.smoothing, 'tikhonov') + strcmp(invert.smoothing, 'cv') + strcmp(invert.smoothing, 'SVD') < 1;                                % if one of these is true, strcmp = 1
    error('ERROR - smoothing must either be ''none'', ''laplacian'' or ''VK'' or ''minimumnorm''. Or ''tikhonov'' or ''cv'' for least squares. Or ''SVD''. What were YOU asking me to do?!');
end


if strcmp(priors.slip_prior, 'boxcar') + strcmp(priors.slip_prior, 'gaussian') + strcmp(priors.slip_prior, 'logarithmic') < 1                           % if one of these is true, strcmp = 1
    error('ERROR - slip_prior must either be ''boxcar'' or ''gaussian'' or ''logarithmic''. Honestly. What are you like?');
end

if strcmp(priors.slip_prior, 'gaussian') == 1                    
    error('ERROR - Unfortunately, slipBERI doesn''t have the capability to draw from a gaussian slip prior. Coming soon!! Maybe.');
end

if strcmp(invert.variable_rake, 'yes') + strcmp(invert.variable_rake, 'no') < 1                                        % if one of these is true, strcmp = 1
    error('ERROR - variable_rake must either be ''yes'' or ''no''. Seriously, what''s your problem?');          
end

if priors.min_slip > priors.max_slip
   error('you''ve got your permitted min_slip and max_slip values the wrong way round') 
end

%if priors.min_dip > priors.max_dip
%   error('you''ve got your permitted min_dip and max_dip values the wrong way round') 
%end

if priors.min_alpha2 > priors.max_alpha2 | priors.min_alpha2 > priors.max_alpha2
   error('you''ve got your permitted min_alpha2 and max_alpha2 values the wrong way round') 
end

if priors.min_rake > priors.max_rake
   error('you''ve got your permitted min_rake and max_rake values the wrong way round') 
end

if priors.min_offset > priors.max_offset
   error('you''ve got your permitted min_offset and max_offset values the wrong way round') 
end

if priors.min_beta > priors.max_beta
   error('you''ve got your permitted min_beta and max_beta values the wrong way round') 
end

if sum(strcmp(data.InSAR_datafile, 'none')) ~= 1 && strcmp(data.InSAR_coordinate_unit, 'utm') + strcmp(data.InSAR_coordinate_unit, 'long/lat') < 1                           % if one of these is true, strcmp = 1
    error('ERROR - InSAR_coordinate_unit must either be ''utm'' or ''long/lat''');
end

if strcmp(testing.testing_mode, 'no')
    if strcmp(data.InSAR_datafile, 'none') + strcmp(data.GPS_datafile_2d, 'none') + strcmp(data.GPS_datafile_3d, 'none') + strcmp(data.atolls_datafile, 'none') == 4  % if one of these is true, strcmp = 1
        error('You don''t have any datasets?!');
    end
end

if strcmp(invert.ensemble_sampling, 'yes') + strcmp(invert.ensemble_sampling, 'no') ~= 1
    error('Whether to use ensemble sampling must be ''yes'' or ''no'' in invert.ensemble_sampling')
end

if strcmp(invert.inversion_type, 'SVD') == 1 && strcmp(invert.smoothing, 'none') ~= 1
    disp('Since you''re doing an SVD inversion, smoothing has been changed to ''none'' and iterations has been changed to ''1''')
    invert.smoothing = 'none';
    invert.iterations = 1;
end

if strcmp(invert.inversion_type, 'least_squares') ==1 && strcmp(invert.smoothing, 'VK') == 1
    error('I''m terribly sorry, slipBERI can''t solve for VK smoothing using a least_squares inversion at this moment in time. Would you care to choose something else?')
end

if strcmp(invert.solve_for_correlation_length, 'yes') ==1 && strcmp(invert.smoothing, 'laplacian') == 1
   disp('')
   disp('Solving for correlation lengths is only for VK smoothing, have changed invert.solve_for_correlation_length from ''yes'' to ''no''');
   invert.solve_for_correlation_length = 'no';
end

if strcmp(data.InSAR_datafile, 'none') == 1
    data.InSAR_coordinate_unit = [];
    invert.solve_for_InSAR_offset = 'no';
end

if strcmp(invert.solve_for_InSAR_offset, 'yes') + strcmp(invert.solve_for_InSAR_ramp, 'yes') == 2
   disp('')
   disp('Solving for InSAR ramp instead of InSAR offset')
   disp('invert.solve_for_InSAR_offset changed to no')
   disp('invert.solve_for_InSAR_ramp left as yes')
   disp('')
   invert.solve_for_InSAR_offset = 'no';
end

if strcmp(data.atolls_datafile, 'none') == 1
    data.atolls_coordinate_unit = [];
end

if strcmp(data.GPS_datafile_2d, 'none') == 1 && strcmp(data.GPS_datafile_3d, 'none') == 1
    data.GPS_coordinate_unit = [];
end

if data.weight_InSAR - data.weight_GPS == 0 && data.weight_InSAR ~= 1
    disp('If you want to equally weight your InSAR and GPS, then data.weight_GPS and data.weight_InSAR both must equal 1');
    disp('(and not any other number.) I''ve changed it so data.weight_GPS=1 and data.weight_InSAR=1');
    m=input('Is that okay? Type y or n and hit enter: ','s');
    if m=='y'
        data.weight_InSAR = 1;
        data.weight_GPS = 1;
        disp('Sorted.');
    elseif m=='n'
        keyboard
    end
end


%% Solving for fault size
% if strcmp(invert.solve_for_fault_size, 'yes') == 1 && strcmp(invert.circular_harmonics, 'no') == 1
%     disp('If solving for fault size, you must choose circular harmonics.')
%     disp('So you need ''invert.circular_harmonics'' to be ''yes''.');
%     disp('I''ve changed this for you')
%     invert.circular_harmonics = 'yes';
% end

if strcmp(invert.solve_for_fault_size, 'yes') == 1 && isempty(priors.min_circharm_phi) == 1 || strcmp(invert.solve_for_fault_size, 'yes') == 1 && isempty(priors.max_circharm_phi) == 1 || strcmp(invert.solve_for_fault_size, 'yes') == 1 && isempty(priors.max_circharm_coeffs) == 1 || strcmp(invert.solve_for_fault_size, 'yes') == 1 && isempty(priors.min_circharm_coeffs) == 1 || strcmp(invert.solve_for_fault_size, 'yes') == 1 && isempty(priors.min_circharm_center) == 1 || strcmp(invert.solve_for_fault_size, 'yes') == 1 && isempty(priors.max_circharm_center) == 1
   error('Circharm priors cannot be empty') 
end


if priors.min_circharm_coeffs > priors.max_circharm_coeffs
   error('you''ve got your permitted min_circharm_coeffs and max_circharm_coeffs values the wrong way round') 
end

if priors.min_circharm_phi > priors.max_circharm_phi
   error('you''ve got your permitted min_circharm_phi and max_circharm_phi values the wrong way round') 
end

if priors.min_circharm_center > priors.max_circharm_center
   error('you''ve got your permitted min_circharm_center and max_circharm_center values the wrong way round') 
end

if strcmp(testing.testing_mode, 'yes') == 1
    data.InSAR_datefile = 'none';
    data.GPS_datefile_2d = 'none';
    data.GPS_datefile_3d = 'none';
end

%% Ensemble sampling
if strcmp(invert.ensemble_sampling, 'no') == 1
   invert.ensemble_move_style = [];
   invert.ensemble_start = [];
end

if strcmp(invert.ensemble_sampling, 'yes') == 1 && strcmp(invert.ensemble_move_style, 'walk') + strcmp(invert.ensemble_move_style, 'stretch') ~=1
   error('You must choose either ''walk'' or ''stretch'' for invert.ensemble_move_style.') 
end


if strcmp(invert.ensemble_sampling, 'yes') == 1 && strcmp(invert.ensemble_start, 'scatter') + strcmp(invert.ensemble_start, 'tight') ~= 1
    error('Select whether you want invert.ensemble_start to be ''tight'' around your initial parameters or ''scatter''')
end

if exist('invert.offset_initial') == 0
    invert.offset_initial = 0;
end





disp('Checks successful.') 