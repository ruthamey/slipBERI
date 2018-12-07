function [] = SSE_inversion_model_visualisation(type_flag,interface_type,VAR_COVAR_flag,smoothness_flag,number_files,save_fig_flag,n_big_patches_keep)
% this is the last version
% 24/09/2013    DB  Make this version compatible with leeds cluster
% 25/09/2013    DB  Add all the option paths
% 25/09/2013    DB  Check the insar plance scaling correction! Adding units
%                   to the histogram plot
% 25/09/2013    DB  Adding the optimum smoothness model and plot results
% 18/11/2014    DB  Add preprocessed datapaths


% save_fig_flag = 0;              % when 1 saving the figures
% type_flag = 'after';           % 'before', 'after' or 'GPS']
% interface_type = 1;             % interface either 1 or 2 
% VAR_COVAR_flag = 'COVAR';         % 'VAR' or 'COVAR'. Does not matter for GPS
% number_files = 250;             % total number of file simulations
% n_big_patches_keep = 1;         % Make the prob distr of the X bigest slipping patches
close all


plot_interface = 1;             % 1 to plot and 0 to supress plotting of the interface
skip_files = 100;                 % skip x first files

old_roughness_method = 0;

checkerboard = 'n'
checkerboard_folder = 'no_noise';
checkerboard_group_number = 3;

if isempty(smoothness_flag)
    plot_slip_magnitudes = 1;
else
     plot_slip_magnitudes = 0;   
end

% file paths
if strcmp(checkerboard,'y')
    
    keyboard
    checkerboard_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/Checkerboard_test' filesep checkerboard_folder filesep 'Interface_' num2str(interface_type) filesep 'group_' num2str(checkerboard_group_number)];
    
    if strcmp(type_flag,'GPS')
        GPS_alone = 1;
        outdata_path = [checkerboard_path filesep 'output' filesep 'GPS_250.mat'];
        outdata_path_files =  [checkerboard_path filesep 'output' filesep 'GPS' ];
        if interface_type ==1
            model_greens_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw/Interface_1/okada_greens_caltec_3stages_high.mat' ];
        else
           model_greens_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw/Interface_2/okada_greens_Mathilde_high.mat' ];
        end
    else
        GPS_alone = 0;
        outdata_path = [checkerboard_path filesep 'output/run_21_Jan_2014/' filesep VAR_COVAR_flag '_250.mat'];
        outdata_path_files =  [checkerboard_path filesep 'output/run_21_Jan_2014/' filesep VAR_COVAR_flag ];
    end

else
    if strcmp(type_flag,'after')
        if strcmp(smoothness_flag,'smooth')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp_smooth/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp_smooth/' VAR_COVAR_flag ];
            GPS_alone = 0;
        elseif strcmp(smoothness_flag,'rough')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rough/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rough/' VAR_COVAR_flag ];
            GPS_alone = 0;
        elseif strcmp(smoothness_flag,'rougher')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rougher/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rougher/' VAR_COVAR_flag ];
            GPS_alone = 0;
        else
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_powerlaw_new_fix_error/Interface_' num2str(interface_type) '/output/temp/' VAR_COVAR_flag ];
            GPS_alone = 0;
        end
    elseif strcmp(type_flag,'before')
        if strcmp(smoothness_flag,'smooth')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_smooth/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_smooth/' VAR_COVAR_flag ];
            GPS_alone = 0;
        elseif strcmp(smoothness_flag,'rough')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rough/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rough/' VAR_COVAR_flag ];
            GPS_alone = 0;
        elseif strcmp(smoothness_flag,'rougher')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rougher/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rougher/' VAR_COVAR_flag ];
            GPS_alone = 0;
        else
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp/' VAR_COVAR_flag '_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp/' VAR_COVAR_flag ];
            GPS_alone = 0;
        end
    elseif  strcmp(type_flag,'GPS')
        if strcmp(smoothness_flag,'smooth')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_smooth/GPS_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_smooth/GPS'  ];
        elseif strcmp(smoothness_flag,'rough')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rough/GPS_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rough/GPS'  ];
        elseif strcmp(smoothness_flag,'rougher')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rougher/GPS_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_rougher/GPS'  ];
        elseif strcmp(smoothness_flag,'roughest')
            outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_roughest/GPS_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp_roughest/GPS'  ];
        else
           outdata_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp/GPS_250.mat'];
            outdata_path_files =  ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_' num2str(interface_type) '/output/temp/GPS'  ];
        end
        if interface_type ==1
            model_greens_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_1/okada_greens_caltec_3stages_high.mat' ];
        else
           model_greens_path = ['/nfs/a1/insar/mexico/envisat/track_255/SSE_model/CASE_tca_no_new_fix_error/Interface_2/okada_greens_Mathilde_high.mat' ];
        end
        GPS_alone = 1;
    end
end

Mexico_borders_path = '/nfs/a1/insar/mexico/envisat/track_255/SSE_model/programs/MEX_adm';
InSAR_border_file = '/nfs/a1/insar/mexico/envisat/track_255/SSE_model/programs/InSAR_box.mat';
rupture_zone_file = '/nfs/a1/insar/mexico/envisat/track_255/SSE_model/programs/rupture_zones/rupture_zones.mat';


outdata_path
% plotting options
marker_size = 75;
rake_hist = 1;
plane_hist = 1;
R_InSAR = 10000;         % radius for the InSAR profile cut

flag_clims = 1;         % if set to 1 use the color limits as specified below.
clims_obs = [-3 6];
clims_est_obs = [-3 6]; 
clims_res = [-3 6];
clims_slip = [0 0.2] ;
clims_slip_std_rake = [0 0.1];
clims_std_rake = [0 20]
% information of the model code
sigma_exclude = 1;          % when 1 sigma was excluded in the simulation run
smoothness_exclude = [];     % when 1 there was no smoothness assumed
USL_flag=1;
NVT_flag = 1;

%% figure properties

fontsize  = 18;
fontsize_GPS = 30;
% GPS vectors
displ_unit_scale = 5;         % [cm] corresonds to legend arrow 
displ_scale = 0.1;    
displ_unit_scale_res = 1;         % [cm] corresonds to legend arrow 
displ_scale_res = 0.5;
displ_color = 'k';
displ_width = 2;



% cross-sectional points:
point2 = {[-99.1 19.75],[-98.8225   19.7217],[ -98.5654   19.6710 ]};
point1 = {[-99.8 16.6],[  -99.5312   16.5486],[  -99.2395   16.5119]};
 
   

% ylimits = [15.5 20];
% xlimits = [-102 -96];
ylimits = [15.5 20];
xlimits = [-102.5 -96.5];
% limnits for the surface deformations
ylimits2 = [16.25 20];
xlimits2 = [-100.75 -98];

%% No changes below:
figprop.xlimits = xlimits;
figprop.ylimits = ylimits ;
figprop.xlimits2 = xlimits2;
figprop.ylimits2 = ylimits2 ;
figprop.clims_slip =clims_slip;
figprop.clims_obs = clims_obs;
figprop.clims_est_obs = clims_est_obs;
figprop.clims_res = clims_res;
figprop.USL_flag = USL_flag;
figprop.NVT_flag = NVT_flag;
figprop.plot_slip_magnitudes=plot_slip_magnitudes;
figprop.fontsize =fontsize;
figprop.flag_clims=flag_clims;
figprop.displ_unit_scale = displ_unit_scale;
figprop.displ_unit_scale_res =displ_unit_scale_res
figprop.displ_scale = displ_scale;
figprop.displ_scale_res = displ_scale_res;
figprop.displ_color = displ_color;
figprop.displ_width = displ_width;
figprop.marker_size = marker_size;
figprop.point1=point1;
figprop.point2=point2;
figprop.R_InSAR=R_InSAR;



[temp1, name, temp2] = fileparts(outdata_path);

clear temp1 temp2
if GPS_alone==0 && sum(name(1:3)=='GPS')==3
    GPS_alone=1;
    fprintf('Seems to be GPS alone. If this is the case set GPS_alone to 1. \n')
    keyboard
end
clear name


if exist(outdata_path_files,'dir')==0
   mkdir(outdata_path_files) 
else
    if save_fig_flag==1
        fprintf('Delete the exisiting files and save new ones \n')
       delete([outdata_path_files filesep '*'])
    end
end

%% Visualisation of the defined interfaces
if plot_interface==1
    subduction_interface_visualisation
end



%% Loading data of the final model
outdata = load(outdata_path);
cal0 = outdata.cal0;
m = outdata.model_interface;
m_opt = outdata.m_opt;
n_patch = outdata.n_patch;
smoothnes_fixed_flag =0;
if isfield(outdata,'with_smoothness')
    if strcmp(outdata.with_smoothness,'y')
        smoothness_exclude =0;
    elseif strcmp(outdata.with_smoothness,'n')
        smoothness_exclude =1;
        ix_big_patches = [62 99 129];
    elseif strcmp(outdata.with_smoothness,'smooth')
        smoothness_exclude =0;
        smoothnes_fixed_flag =1;
    elseif strcmp(outdata.with_smoothness,'rough')
        smoothness_exclude =0;
        smoothnes_fixed_flag =1;

        
    end
else
    smoothness_exclude =0;
end
if isfield(outdata,'transitions')==1
    transitions = outdata.transitions;
else
    fprintf('interface transitions are not plotted! \n')
    transitions=[];
end

% GPS observations
obs_GPS_opt= outdata.obs_GPS_opt;
obs_GPS = outdata.obs_GPS;
% The covariance matrices used in the estimation
Q_GPS = outdata.Q_GPS;
% GPS longitude and latitude
ll_GPS = outdata.ll_GPS;
% number of GPS observations
n_GPS = size(ll_GPS,1);
% Decomposing GPS modelled observations vector into multiple components 
obs_GPS_opt_E = obs_GPS_opt(1:n_GPS);
obs_GPS_opt_N = obs_GPS_opt(n_GPS+1:2*n_GPS);
obs_GPS_opt_U = obs_GPS_opt(2*n_GPS+1:end);
% Decomposing GPS observations vector into multiple components 
obs_GPS_E = obs_GPS(1:n_GPS);
obs_GPS_N = obs_GPS(n_GPS+1:2*n_GPS);
obs_GPS_U = obs_GPS(2*n_GPS+1:end);

% InSAR observations
if GPS_alone~=1
    % GPS and InSAR observations
    obs_SAR_opt= outdata.obs_SAR_opt;
    obs_SAR = outdata.obs_SAR;
    % The covariance matrices used in the estimation
    Q_SAR = outdata.Q_SAR;
    % Estimated InSAR plane
    obs_SAR_plane = outdata.obs_SAR_plane;
    % GPS and SAR longitude and latitude
    ll_SAR = outdata.ll_SAR;
    % InSAR observations and modelled observations corrected for plane
    obs_SAR_min_plane = obs_SAR - obs_SAR_plane;
    obs_SAR_opt_min_plane = obs_SAR_opt - obs_SAR_plane;
end

% Getting the slip from the outfile
rake=m_opt(1:n_patch);
rake_slip = m_opt(n_patch+1:2*n_patch);
strike_slip =m_opt(n_patch+1:2*n_patch).*cos(rake*pi/180);
dip_slip = m_opt(n_patch+1:2*n_patch).*sin(rake*pi/180);


% saving the data for Sylvain 
temp_path = fileparts(outdata_path);
temp_pathfile = [temp_path filesep 'model.mat'];
save(temp_pathfile,'rake_slip','rake')

keyboard
if smoothness_exclude~=1
    %% Plotting the profile of slip
    [temp1,temp2,interface_center_vertices] = drawmodel(m,'origin',cal0);
    interface.depth = interface_center_vertices(3,:)'; 
    interface.lonlat = interface_center_vertices(1:2,:)';
    clear interface_center_vertices temp1 temp2
    interface.rake_slip =rake_slip; 
    interface.strike_slip =strike_slip; 
    interface.dip_slip =dip_slip; 

    point1_interface = [-99.8570 16.8220];
    point2_interface = [-99.1 19.75];
    point3_interface = [ -100.1820   15.5365];
    x_limits = [-100 325];


    [interface_profile.rake_slip interface_profile.XY]  = profile_interpolation(interface.rake_slip,interface.lonlat,point1_interface,point2_interface);
    [interface_profile.strike_slip interface_profile.XY]  = profile_interpolation(interface.strike_slip,interface.lonlat,point1_interface,point2_interface);
    [interface_profile.dip_slip interface_profile.XY]  = profile_interpolation(interface.dip_slip,interface.lonlat,point1_interface,point2_interface);


    h1 = figure('position',[680   697   814   381]);
    set(gca,'position',[    0.1300    0.1365    0.7750    0.7885])
    plot(interface_profile.XY(:,1),interface_profile.dip_slip,':','color',[0.65 0.65 0.65],'linewidth',4)
    hold on
    plot(interface_profile.XY(:,1),interface_profile.strike_slip,'b--','linewidth',4)
    plot(interface_profile.XY(:,1),interface_profile.rake_slip,'r.-','linewidth',4)
    hold off
    legend('dip slip','strike slip','rake slip',0)
    legend boxoff
    xlim(x_limits)
    ylabel('slip [m]','fontsize',fontsize)
    xlabel('Horizontal distance from ACAP [km]','fontsize',fontsize)
    set(gca,'fontsize',fontsize)
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'slip_profile.eps'])
    end
    clear h1


    %% Visualisations

    data.dip_slip = dip_slip;
    data.strike_slip = strike_slip;
    data.rake_slip = rake_slip;
    data.m =  outdata.model_interface;
    data.ll_GPS = ll_GPS;
    data.cal0 = cal0;
    data.transitions = transitions;
    data.GPS_alone=GPS_alone;
    data.obs_GPS_E = obs_GPS_E;
    data.obs_GPS_N = obs_GPS_N;
    data.obs_GPS_U = obs_GPS_U;
    data.obs_GPS_opt_E = obs_GPS_opt_E;
    data.obs_GPS_opt_N = obs_GPS_opt_N;
    data.obs_GPS_opt_U = obs_GPS_opt_U;
    data.A_design =  outdata.A_design;
    data.outdata_path_files = outdata_path_files;
    data.GPS_LOS_struct = outdata.obs_GPS_LOS;
    data.Q_GPS = Q_GPS;
    data.checkerboard = checkerboard;
    if GPS_alone~=1
        data.ll_SAR = ll_SAR;
        data.obs_SAR = obs_SAR;
        data.obs_SAR_opt = obs_SAR_opt;
        data.obs_SAR_plane = obs_SAR_plane;
        data.Q_SAR = Q_SAR;
    else
        % give the strike and dip greens parameters of the InSAR such the model
        % can be shown.
        temp_greens = load(model_greens_path);
        data.Gss_SAR = temp_greens.Gss_SAR;
        data.Gds_SAR = temp_greens.Gds_SAR;
        data.ll_SAR = temp_greens.ll_SAR;
    end
    plotting_model(data,figprop,'MaxLikelihood')
end



%% Loading the data from all the simulations in between

if number_files-skip_files>=1
    % identifying the 5 biggest slipping patches
    rake_slip_sorted  = sort(abs(rake_slip));
    if smoothness_exclude==1
        for kk=1:n_big_patches_keep
            ix_big_patches(kk) = find( rake_slip_sorted(end-kk+1) == abs(rake_slip));
        end
        clear kk
    else
        for kk=1:n_big_patches_keep
            ix_big_patches(kk) = find( rake_slip_sorted(end-kk+1) == abs(rake_slip));
        end
        clear kk
    end

    % optimum model parameters
    rake_slip_opt = rake_slip(ix_big_patches);
    rake_opt = rake(ix_big_patches);
    strike_slip_opt = strike_slip(ix_big_patches);
    dip_slip_opt = dip_slip(ix_big_patches);
    
    % cope with two codes for the model inversion
    % following coeffients are for the model and do not vary for the
    % biggest slip patches!
    if sigma_exclude~=1
        if smoothness_exclude==1        % no smoothness
            if GPS_alone~=1
                plane_opt = m_opt(end-3:end-1);
            end
            hyper_opt = m_opt(end);
        else
            if GPS_alone~=1
                plane_opt = m_opt(end-4:end-2);
            end
            hyper_opt = m_opt(end-1:end);
        end
    else
        if smoothness_exclude==1        % no smoothness
            if GPS_alone~=1
                plane_opt = m_opt(end-2:end);
            end
            hyper_opt = []; 
        else
            if GPS_alone~=1
                plane_opt = m_opt(end-3:end-1);
            end
            hyper_opt = m_opt(end);  
        end
    end

    % Keep only data of the 5 biggest slipping patches, the plane and the hyperparameters
    % Keep all the rake slip in order to make the standard deviation plot
    for k=1:number_files-skip_files
        outdata_temp = load([outdata_path_files '_' num2str(skip_files+k) '.mat']);
       
        if k==1
           fprintf('Loading files \n')
           n_runs_file =  size(outdata_temp.m_keep,2);
           % initialisation of variables
           dip_slip_keep = zeros([length(ix_big_patches) (number_files-skip_files)*n_runs_file]);
           strike_slip_keep = zeros([length(ix_big_patches) (number_files-skip_files)*n_runs_file]);
           if GPS_alone~=1
                plane_keep = zeros([3 (number_files-skip_files)*n_runs_file]);        % [a x + by  +c], matrix size = [3_coeff n_sim ]
           end
           % cope with two codes for the model inversion
           if sigma_exclude~=1
               if smoothness_exclude==1
                   % eclude smoothness
                  hyper_keep = zeros([1 (number_files-skip_files)*n_runs_file]);        % [sigma]
               else
                   hyper_keep = zeros([2 (number_files-skip_files)*n_runs_file]);        % [smoothness ; sigma]
               end
           else
               if smoothness_exclude==1
                    % exlude smoothness
                    hyper_keep = [];
               else
                   hyper_keep = zeros([1 (number_files-skip_files)*n_runs_file]);        % [smoothness]
               end
           end
           
           % Keep of all the rake slip
           rake_slip_all = single(zeros([n_patch (number_files-skip_files)*n_runs_file]));
           rake_all = single(zeros([n_patch (number_files-skip_files)*n_runs_file]));

        end
        if k == floor(k/10)*10
            fprintf([num2str(k) ' of ' num2str(number_files-skip_files) '\n'])
        end
        % m is a vector which contains [rake for each patch, slip for each patch, plane 3 coeff ,hyperparameter alpha, (hyperparameter sigma)]
        m_keep_temp = outdata_temp.m_keep;


        % Probability keep
        P_keep_all(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = outdata_temp.P_keep;
        P_post_keep_all(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = outdata_temp.P_post_keep;
        % cope with two codes for the model inversion
        if sigma_exclude~=1
            if smoothness_exclude==1  % no smoothness
                hyper_keep(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = m_keep_temp(end,:);
                if GPS_alone~=1
                    % extracting the plane coefficients
                    plane_keep(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = m_keep_temp(end-3:end-1,:);                  
                end
            else
                % extracting the hyperparameter and the smoothness
                hyper_keep(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = m_keep_temp(end-1:end,:);
                if GPS_alone~=1
                    % extracting the plane coefficients
                    plane_keep(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = m_keep_temp(end-4:end-2,:);                  
                end
            end
        else
            if smoothness_exclude==1
                hyper_keep=[];
                if GPS_alone~=1
                    % extracting the plane coefficients
                    plane_keep(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = m_keep_temp(end-2:end,:);
                end
            else
                % extracting the smoothness            
                hyper_keep(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = m_keep_temp(end,:);
                if GPS_alone~=1
                    % extracting the plane coefficients
                    plane_keep(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = m_keep_temp(end-3:end-1,:);
                end
            end
        end
        % extracting the dip/strike slip from the rake and slip 
        rake_temp=m_keep_temp(1:n_patch,:);
        rake_keep(1:n_big_patches_keep,n_runs_file*(k-1)+1:n_runs_file*(k)) = rake_temp(ix_big_patches,:);
        rake_keep_all(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = rake_temp;
        
        rake_slip_temp=m_keep_temp(n_patch+1:n_patch*2,:);
        rake_slip_keep (1:n_big_patches_keep,n_runs_file*(k-1)+1:n_runs_file*(k))= rake_slip_temp(ix_big_patches,:);
        rake_slip_keep_all (:,n_runs_file*(k-1)+1:n_runs_file*(k))= rake_slip_temp;

        clear m_keep_temp
        strike_slip_temp = rake_slip_temp.*sin(rake_temp*pi/180);          % strike slip
        dip_slip_temp= rake_slip_temp.*cos(rake_temp*pi/180);              % dip slip
        
        % extracting all the rake slip
        rake_slip_all(:,n_runs_file*(k-1)+1:n_runs_file*(k))=single(rake_slip_temp);
        rake_all(:,n_runs_file*(k-1)+1:n_runs_file*(k))=single(rake_temp);
        clear rake_slip_temp rake_temp
        % keeping only 5 biggest slipping patches
        dip_slip_keep(1:n_big_patches_keep,n_runs_file*(k-1)+1:n_runs_file*(k)) = dip_slip_temp(ix_big_patches,:);
        dip_slip_keep_all(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = dip_slip_temp;

        clear dip_slip_temp
        strike_slip_keep(1:n_big_patches_keep,n_runs_file*(k-1)+1:n_runs_file*(k)) = strike_slip_temp(ix_big_patches,:);
        strike_slip_keep_all(:,n_runs_file*(k-1)+1:n_runs_file*(k)) = strike_slip_temp;

        clear strike_slip_temp
    end
    
    if GPS_alone~=1
        % rescalling the plane coefficients
        plane_keep = plane_keep./100000.*100;        % plane units are now in cm
        plane_opt = plane_opt./100000.*100;        % plane units are now in cm
    end
    
    % changing to the right scale of the hyperparameters
    hyper_keep = 10.^(hyper_keep);
    hyper_opt = 10.^(hyper_opt);


    % resetting the plotting matrix options:
    if GPS_alone==1
        plane_hist = 0;
    end
    
    
    if smoothnes_fixed_flag==1
        if round((mean(hyper_keep)-hyper_opt)*1000)/1000~=0
            fprintf('something seems wrong with the smoothness \n')
            keyboard
        end
        hyper_keep =[];
        hyper_opt = [];
    end
    
    
    
    %% Plotting the masked maximum likelihood solution    
    rake_std = std(rake_slip_keep_all');
    rake_slip_lower_bound = data.rake_slip-2*rake_std';
    ix=find(rake_slip_lower_bound<0);

    data_masked= data;
    data_masked.rake_slip(ix) =0; 
    data_masked.strike_slip(ix) =0; 
    data_masked.dip_slip(ix) =0; 

    plotting_model(data_masked,figprop,'MaxLikelihood_masked')
    
 
    
    %% Plotting the mean of the distribution
    if number_files-skip_files>=1

            figprop_mean = figprop;
            figprop_mean.colorbar_unit='m';

            % computation of the mean of the rake and the rake slip 
            rake_slip_mean = mean(rake_slip_all');
            % computation of the mean of the rake and the rake slip 
            rake_mean = mean(rake_all');

            strike_slip_mean = rake_slip_mean.*sin(rake_mean*pi/180);          % strike slip
            dip_slip_mean= rake_slip_mean.*cos(rake_mean*pi/180);              % dip slip
            plot_data = rake_slip_mean;

            % setting up the data matrix
            data_mean.plot_data_slipcolor = 1;
            data_mean.plot_data = plot_data;
            data_mean.m =  m;
            data_mean.m(8,:) = strike_slip_mean;
            data_mean.m(9,:) = dip_slip_mean;
            data_mean.ll_GPS = ll_GPS;
            data_mean.cal0 = cal0;  
            data_mean.transitions = transitions;
            data_mean.GPS_alone=GPS_alone;
            data_mean.outdata_path_files = outdata_path_files;

            if GPS_alone~=1
                data_mean.ll_SAR = ll_SAR;
            else
                temp_greens = load(model_greens_path);
                data_mean.ll_SAR = temp_greens.ll_SAR;
            end

            % plotting the results
            plotting_model(data_mean,figprop_mean,'Mean slip solution') 

    end
    
    
    if smoothness_exclude==1 
        % Replace the maximum likelihood by the mean
        dip_slip_opt = dip_slip_mean(ix_big_patches)';
        strike_slip_opt = strike_slip_mean(ix_big_patches)';
        rake_opt = rake_mean(ix_big_patches)';
        rake_slip_opt = rake_slip_mean(ix_big_patches)';
    end
    


    % [H,AX,BigAx,P,PAx] = PLOTMATRIX(...) returns a matrix of handles
    % to the objects created in H, a matrix of handles to the individual
    % subaxes in AX, a handle to big (invisible) axes that frame the
    % subaxes in BigAx, a matrix of handles for the histogram plots in
    % P, and a matrix of handles for invisible axes that control the
    % histogram axes scales in PAx.  BigAx is left as the CurrentAxes so
    % that a subsequent TITLE, XLABEL, or YLABEL will be centered with
    % respect to the matrix of axes.

    % For the scatter plots
    for kk=1:n_big_patches_keep
        if rake_hist~=1
            if plane_hist==1
                data = [dip_slip_keep(kk,:)', strike_slip_keep(kk,:)', hyper_keep', plane_keep'];
                data_opt = [dip_slip_opt(kk,:)', strike_slip_opt(kk,:)', hyper_opt', plane_opt'];
                if sigma_exclude~=1
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1
                        strings{1} = 'Dip';                 strings_units{1} = 'Dip [m]';
                        strings{2} = 'Strike';              strings_units{2} = 'Strike [m]';
                        strings{3} = '\sigma';              strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{4} = 'a';                   strings_units{4} = 'a [cm/km]';
                        strings{5} = 'b';                   strings_units{5} = 'b [cm/km]';
                        strings{6} = 'c';                   strings_units{6} = 'c [cm]';
                        strings{7} = 'Frequency';
                    else
                        strings{1} = 'Dip';                 strings_units{1} = 'Dip [m]';
                        strings{2} = 'Strike';              strings_units{2} = 'Strike [m]';
                        strings{4} = '\sigma';              strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{3} = '\alpha^2';              strings_units{4} = 'Smoothness \alpha^2    ';
                        strings{5} = 'a';                   strings_units{5} = 'a [cm/km]';
                        strings{6} = 'b';                   strings_units{6} = 'b [cm/km]';
                        strings{7} = 'c';                   strings_units{7} = 'c [cm]';
                        strings{8} = 'Frequency';
                        end
                else
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1
                            strings{1} = 'Dip';                 strings_units{1} = 'Dip [m]';
                            strings{2} = 'Strike';              strings_units{2} = 'Strike [m]';
                            strings{3} = 'a';                   strings_units{3} = 'a [cm/km]';
                            strings{4} = 'b';                   strings_units{4} = 'b [cm/km]';
                            strings{5} = 'c';                   strings_units{5} = 'c [cm]';
                            strings{6} = 'Frequency';
                    else
                        strings{1} = 'Dip';                 strings_units{1} = 'Dip [m]';
                            strings{2} = 'Strike';              strings_units{2} = 'Strike [m]';
                            strings{3} = '\alpha^2';              strings_units{3} = 'Smoothness \alpha^2    ';
                            strings{4} = 'a';                   strings_units{4} = 'a [cm/km]';
                            strings{5} = 'b';                   strings_units{5} = 'b [cm/km]';
                            strings{6} = 'c';                   strings_units{6} = 'c [cm]';
                            strings{7} = 'Frequency';
                    end
                end
            else
                data = [dip_slip_keep(kk,:)', strike_slip_keep(kk,:)', hyper_keep'];
                data_opt = [dip_slip_opt(kk,:)', strike_slip_opt(kk,:)', hyper_opt'];
                if sigma_exclude~=1
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1
                        strings{1} = 'Dip';                  strings_units{1} = 'Dip [m]';
                        strings{2} = 'Strike';               strings_units{2} = 'Strike [m]';
                        strings{3} = '\sigma';               strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{4} = 'Frequency';
                    else
                        strings{1} = 'Dip';                  strings_units{1} = 'Dip [m]';
                        strings{2} = 'Strike';               strings_units{2} = 'Strike [m]';
                        strings{4} = '\sigma';               strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{3} = '\alpha^2';               strings_units{4} = 'Smoothness \alpha^2    ';
                        strings{5} = 'Frequency';
                    end
                else
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1
                        strings{1} = 'Dip';                  strings_units{1} = 'Dip [m]';
                        strings{2} = 'Strike';               strings_units{2} = 'Strike [m]';
                        strings{3} = 'Frequency'; 
                    else
                        strings{1} = 'Dip';                  strings_units{1} = 'Dip [m]';
                        strings{2} = 'Strike';               strings_units{2} = 'Strike [m]';
                        strings{3} = '\alpha^2';               strings_units{3} = 'Smoothness \alpha^2   ';
                        strings{4} = 'Frequency'; 
                    end
                end
            end
        else
            if plane_hist==1
                data = [rake_slip_keep(kk,:)', rake_keep(kk,:)', hyper_keep', plane_keep'];
                data_opt = [rake_slip_opt(kk,:)', rake_opt(kk,:)', hyper_opt', plane_opt'];
                if sigma_exclude~=1
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [m]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{3} = '\sigma';               strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{4} = 'a';                    strings_units{4} = 'a [cm/km]';
                        strings{5} = 'b';                    strings_units{5} = 'b [cm/km]';
                        strings{6} = 'c';                    strings_units{6} = 'c [cm]';
                        strings{7} = 'Frequency';
                    else
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [m]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{4} = '\sigma';               strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{3} = '\alpha^2';               strings_units{4} = 'Smoothness \alpha^2    ';
                        strings{5} = 'a';                    strings_units{5} = 'a [cm/km]';
                        strings{6} = 'b';                    strings_units{6} = 'b [cm/km]';
                        strings{7} = 'c';                    strings_units{7} = 'c [cm]';
                        strings{8} = 'Frequency';
                    end
                else
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1 &&  strcmp(checkerboard,'y') 
                        % putting in the log domain
                        data_opt(:,1) = log10(data_opt(:,1));
                        data(:,1) = log10(data(:,1));  
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [log(m)]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{3} = 'a';                    strings_units{3} = 'a [cm/km]';
                        strings{4} = 'b';                    strings_units{4} = 'b [cm/km]';
                        strings{5} = 'c';                    strings_units{5} = 'c [cm]';
                        strings{6} = 'Frequency'; 
                    elseif smoothness_exclude==1 ||  smoothnes_fixed_flag==1 &&  strcmp(checkerboard,'n') 
                        % putting in the log domain
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [m]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{3} = 'a';                    strings_units{3} = 'a [cm/km]';
                        strings{4} = 'b';                    strings_units{4} = 'b [cm/km]';
                        strings{5} = 'c';                    strings_units{5} = 'c [cm]';
                        strings{6} = 'Frequency'; 
                    else
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [m]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{3} = '\alpha^2';               strings_units{3} = 'Smoothness \alpha^2    ';
                        strings{4} = 'a';                    strings_units{4} = 'a [cm/km]';
                        strings{5} = 'b';                    strings_units{5} = 'b [cm/km]';
                        strings{6} = 'c';                    strings_units{6} = 'c [cm]';
                        strings{7} = 'Frequency';
                    end
                end
            else
                data = [rake_slip_keep(kk,:)', rake_keep(kk,:)', hyper_keep'];
                data_opt = [rake_slip_opt(kk,:)', rake_opt(kk,:)', hyper_opt'];
                
                if sigma_exclude~=1
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1
                        % putting in the log domain
                        data_opt(:,1) = log10(data_opt(:,1));
                        data(:,1) = log10(data(:,1));  
                        strings{1} = 'Slip';                            strings_units{1} = 'Slip [log(m)]';
                        strings{2} = 'Rake';                            strings_units{2} = 'Rake [\circ]';
                        strings{3} = '\sigma';                          strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{4} = 'Frequency';
                    else                    
                        strings{1} = 'Slip';                            strings_units{1} = 'Slip [m]';
                        strings{2} = 'Rake';                            strings_units{2} = 'Rake [\circ]';
                        strings{4} = '\sigma';                          strings_units{3} = 'Hyper parameter \sigma [-]';
                        strings{3} = '\alpha^2';                          strings_units{4} = 'Smoothness \alpha^2    ';
                        strings{5} = 'Frequency';
                    end
                else
                    if smoothness_exclude==1 ||  smoothnes_fixed_flag==1  && strcmp(checkerboard,'y')
                        % putting in the log domain
                        data_opt(:,1) = log10(data_opt(:,1));
                        data(:,1) = log10(data(:,1));  
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [log(m)]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{3} = 'Frequency';  
                    elseif smoothness_exclude==1 ||  smoothnes_fixed_flag==1  && strcmp(checkerboard,'n')
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [m]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{3} = 'Frequency'; 
                    else
                        strings{1} = 'Slip';                 strings_units{1} = 'Slip [m]';
                        strings{2} = 'Rake';                 strings_units{2} = 'Rake [\circ]';
                        strings{3} = '\alpha^2';              strings_units{3} = 'Smoothness \alpha^2     ';
                        strings{4} = 'Frequency'; 
                    end
                end
            end
        end
   
        
        if  smoothness_exclude~=1 &&  smoothnes_fixed_flag~=1
            %%% Computing the roughest and maximum smoothness solution
            % roughest
            ix_total_roughness = find(max(hyper_keep(end,:))==hyper_keep(end,:));
            ix_total_roughness = ix_total_roughness(end);
            ix_file_roughness = ceil(ix_total_roughness./n_runs_file)+skip_files;
            % maximum smoothness
            [nelements,xcenters] = hist(hyper_keep(end,:),20);
            ix = find(nelements ==max(nelements));
            binsize = diff(xcenters(1:2));
            ix_total_max_smoothness = find(xcenters(ix)-binsize/2<hyper_keep & hyper_keep<xcenters(ix)+binsize/2);      % the file within all the loaded model runs
            ix_total_max_smoothness = ix_total_max_smoothness(end);                                                     % file which has the highest likelihood
            ix_file_max_smoothness = ceil(ix_total_max_smoothness./n_runs_file)+skip_files;                             % finding the file which contains the model run
            % the solution data
            data_roughness = [data(ix_total_roughness,:)];
            data_max_smoothness = [data(ix_total_max_smoothness,:)];
        end
                
        % plotting the model inversion distributions
        % data = [small slip patch ,....., big slip patch , smoothness , sigma]
        h1= figure('position',[ 680         330        1108         748],'name',['Prob distribution ' num2str(kk) ' biggest slipping patch'])    ;
        set(h1,'renderer','painters')
        if plane_hist==1
            if GPS_alone==1
                set(h1,'Position',[208   371   784   574 ]);
            else
                set(h1,'Position',[216           5        1323         908 ]);
            end
        end
        ix_range = [1:500:size(data,1)];
        temp = unique(ceil(rand(size(data,1),1)*size(data,1)),'stable');
        ix_range_2 = temp(1:ceil(size(data,1)/500));
        ix_range_3 = 1:size(data,1);
        
        if  smoothness_exclude==1 ||  smoothnes_fixed_flag==1 || old_roughness_method~=1
           [H,AX,BigAx,P,PAx,cc_colorbar] = plotmatrix_lower_david(data(ix_range_3,:),'plot_color',data_opt);
        else
           [H,AX,BigAx,P,PAx,cc_colorbar] = plotmatrix_lower_david(data(ix_range_3,:),'plot_color',data_opt,data_max_smoothness,data_roughness);
        end
        
        
        for k=1:length(strings)-1
	    if GPS_alone==1
		fontsize_plot = fontsize_GPS;
	    else
		fontsize_plot = fontsize;
            end
	    set(PAx(k),'fontsize',fontsize_plot-3)
            set(AX(k),'fontsize',fontsize_plot-3)

            set(get(PAx(k),'xlabel'),'string',strings_units{k},'fontsize',fontsize_plot-1)
            set(get(AX(k,1),'ylabel'),'string',strings{k+1},'fontsize',fontsize_plot-1)
            if k==1
                set(get(PAx(k),'ylabel'),'string',strings{end},'fontsize',fontsize_plot-1)
            end
        end
        set(cc_colorbar,'fontsize',fontsize_plot-1)
        title(cc_colorbar,'Frequency','fontsize',fontsize_plot-1)
        set(gcf,'PaperPositionMode','auto')

        
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'prob_dist_' num2str(kk) '_biggest_patch.eps'])
        clear h1
%         if GPS_alone==1
%             keyboard
%         end
    end
end



%% plotting the rake slip standard deviation
if number_files-skip_files>=1


        figprop_rake_std = figprop;
        figprop_rake_std.clims_slip=clims_slip_std_rake;
        figprop_rake_std.colorbar_unit='m';
        
        % computation of the std of the rake slip 
        if smoothness_exclude==1
            
            std_slip_log = std(log10(rake_slip_all'));
            plot_data = std_slip_log;        
            figprop_rake_std.clims_slip=[];
            figprop_rake_std.colorbar_unit='log(m)';


        else
            plot_data = std(rake_slip_all');
        end
        
        % setting up the data matrix
        data_rake_std.plot_data_slipcolor = 1;
        data_rake_std.plot_data = plot_data;
        data_rake_std.m =  outdata_temp.model_interface;
        data_rake_std.ll_GPS = ll_GPS;
        data_rake_std.cal0 = cal0;  
        data_rake_std.transitions = transitions;
        data_rake_std.GPS_alone=GPS_alone;
        data_rake_std.outdata_path_files = outdata_path_files;
     
        if GPS_alone~=1
            data_rake_std.ll_SAR = ll_SAR;
        else
            temp_greens = load(model_greens_path);
            data_rake_std.ll_SAR = temp_greens.ll_SAR;
        end

        % plotting the standard deviations
        plotting_model(data_rake_std,figprop_rake_std,'Slip_std')
        
        
        
        %% standard deviation for the rake
        figprop_rake = figprop;
        figprop_rake.clims_slip=clims_std_rake;
        figprop_rake.colorbar_unit='\circ';
        
        % computation of the std of the rake slip 
        if smoothness_exclude==1
            
            std_rake_log = std(log10(rake_all'));
            plot_data = std_rake_log;        
            figprop_rake.clims_slip=[];
            figprop_rake.colorbar_unit='log(\circ)';


        else
            plot_data = std(rake_all');
        end
        
        % setting up the data matrix
        data_rake_std.plot_data_slipcolor = 1;
        data_rake_std.plot_data = plot_data;
        data_rake_std.m =  outdata_temp.model_interface;
        data_rake_std.ll_GPS = ll_GPS;
        data_rake_std.cal0 = cal0;  
        data_rake_std.transitions = transitions;
        data_rake_std.GPS_alone=GPS_alone;
        data_rake_std.outdata_path_files = outdata_path_files;
     
        if GPS_alone~=1
            data_rake_std.ll_SAR = ll_SAR;
        else
            temp_greens = load(model_greens_path);
            data_rake_std.ll_SAR = temp_greens.ll_SAR;
        end

        % plotting the standard deviations
        plotting_model(data_rake_std,figprop_rake,'Rake_std')
        
       
        %% plotting a masked mean plot in case no smoothness is applied
        if smoothness_exclude==1
            slip_log_range = [min(min(log10(rake_slip_all'))) max(max(log10(rake_slip_all')))];
            uniform_slip_std_log =  sqrt((slip_log_range(2)-slip_log_range(1))^2/12);
            std_threshold = 0.5*uniform_slip_std_log;
            fprintf(['Log slip trhreshold of ' num2str(std_threshold) ' on the standard deviation \n'])
            % finding those that have close to a normal distribution
            ix = find(std_slip_log>std_threshold);
            
            % masking out those which are not constrained
            rake_slip_mean_masked = rake_slip_mean;
            strike_slip_mean_masked = strike_slip_mean;
            dip_slip_mean_masked = dip_slip_mean;
            rake_slip_mean_masked(ix) = NaN;
            strike_slip_mean_masked(ix) =NaN;
            dip_slip_mean_masked(ix) = NaN;
            
            % plotting the results
            
            figprop_mean_masked = figprop;
            figprop_mean_masked.colorbar_unit='m';


            plot_data = rake_slip_mean_masked;

            % setting up the data matrix
            data_mean_masked.plot_data_slipcolor = 1;
            data_mean_masked.plot_data = plot_data;
            data_mean_masked.m =  m;
            data_mean_masked.m(8,:) = strike_slip_mean_masked;
            data_mean_masked.m(9,:) = dip_slip_mean_masked;
            data_mean_masked.ll_GPS = ll_GPS;
            data_mean_masked.cal0 = cal0;  
            data_mean_masked.transitions = transitions;
            data_mean_masked.GPS_alone=GPS_alone;
            data_mean_masked.outdata_path_files = outdata_path_files;

            if GPS_alone~=1
                data_mean_masked.ll_SAR = ll_SAR;
            else
                temp_greens = load(model_greens_path);
                data_mean_masked.ll_SAR = temp_greens.ll_SAR;
            end

            % plotting the results
             % plotting_model(data_mean_masked,figprop_mean_masked,'Mean slip solution masked') 


        end

end




if old_roughness_method==1

    %% Search for the case with min smoothness (roughest solution)
    if number_files-skip_files>=1  && smoothness_exclude~=1   % latter is to check if smoothness was used in the modelling




        % the number of the run within the selected file
        ix_run = ix_total_roughness-(ix_file_roughness-skip_files-1)*n_runs_file;

        % loading the datafile
        outdata_temp = load([outdata_path_files '_' num2str(ix_file_roughness) '.mat']);
        m_opt_max_roughness = outdata_temp.m_keep(:,ix_run);

        % design matrix including greens coefficients and plane
        if GPS_alone~=1
            A_plane = outdata_temp.A_design(:,end-2:end);
            n_sar = length(obs_SAR);
            A_plane_SAR = A_plane(1:n_sar,:);

            if sigma_exclude~=1
                xslope_plane=m_opt_max_roughness(end-4);
                yslope_plane=m_opt_max_roughness(end-3);
                offset_plane=m_opt_max_roughness(end-2);
            else
                xslope_plane=m_opt_max_roughness(end-3);
                yslope_plane=m_opt_max_roughness(end-2);
                offset_plane=m_opt_max_roughness(end-1);
            end

            % Computing the fitted plane to save as output
            obs_SAR_plane_max_roughness = A_plane_SAR*[xslope_plane ;yslope_plane; offset_plane];
        end

        % convert rake slip to dip and strike slip and then ivert using greens
        % to the slip.
        rake_max_roughness=m_opt_max_roughness(1:n_patch);
        rake_slip_max_roughness = m_opt_max_roughness(n_patch+1:2*n_patch);
        strike_slip_max_roughness =m_opt_max_roughness(n_patch+1:2*n_patch).*cos(rake_max_roughness*pi/180);
        dip_slip_max_roughness = m_opt_max_roughness(n_patch+1:2*n_patch).*sin(rake_max_roughness*pi/180);
        m_slip_max_roughness(1:n_patch)=strike_slip_max_roughness;                                    % strike slip
        m_slip_max_roughness(n_patch+1:2*n_patch)=dip_slip_max_roughness;                             % dip slip
        if GPS_alone~=1
            m_slip_max_roughness(2*n_patch+1:2*n_patch+3)=m_opt_max_roughness(2*n_patch+1:2*n_patch+3);   % plane
        end
        clear rake_max_roughness 


        U = outdata_temp.A_design*(m_slip_max_roughness)';
        % splitting surface deformations in a GPS and InSAR component
        if GPS_alone==1
            obs_GPS_opt_max_roughness = U;
        else
            obs_SAR_opt_max_roughness = U(1:n_sar);    
            obs_GPS_opt_max_roughness = U(n_sar+1:end);
        end

        clear U

        % Decomposing GPS modelled observations vector into multiple components 
        obs_GPS_opt_max_roughness_E = obs_GPS_opt_max_roughness(1:n_GPS);
        obs_GPS_opt_max_roughness_N = obs_GPS_opt_max_roughness(n_GPS+1:2*n_GPS);
        obs_GPS_opt_max_roughness_U = obs_GPS_opt_max_roughness(2*n_GPS+1:end);

        % give the strike and dip greens parameters of the InSAR such the model
        % can be shown.
        data_smooth_max_roughness.dip_slip = dip_slip_max_roughness;
        data_smooth_max_roughness.strike_slip = strike_slip_max_roughness;
        data_smooth_max_roughness.rake_slip= rake_slip_max_roughness;
        data_smooth_max_roughness.m =  outdata_temp.model_interface;
        data_smooth_max_roughness.ll_GPS = ll_GPS;
        data_smooth_max_roughness.cal0 = cal0;  
        data_smooth_max_roughness.transitions = transitions;
        data_smooth_max_roughness.GPS_alone=GPS_alone;
        data_smooth_max_roughness.obs_GPS_E = obs_GPS_E;
        data_smooth_max_roughness.obs_GPS_N = obs_GPS_N;
        data_smooth_max_roughness.obs_GPS_U = obs_GPS_U;
        data_smooth_max_roughness.obs_GPS_opt_E = obs_GPS_opt_max_roughness_E;
        data_smooth_max_roughness.obs_GPS_opt_N = obs_GPS_opt_max_roughness_N;
        data_smooth_max_roughness.obs_GPS_opt_U = obs_GPS_opt_max_roughness_U;
        data_smooth_max_roughness.A_design =  outdata_temp.A_design;
        data_smooth_max_roughness.Q_GPS = Q_GPS;
        data_smooth_max_roughness.outdata_path_files = outdata_path_files;
        data_smooth_max_roughness.GPS_LOS_struct = outdata.obs_GPS_LOS;
        data_smooth_max_roughness.checkerboard = checkerboard;

        if GPS_alone~=1
            data_smooth_max_roughness.ll_SAR = ll_SAR;
            data_smooth_max_roughness.obs_SAR = obs_SAR;
            data_smooth_max_roughness.obs_SAR_opt = obs_SAR_opt_max_roughness;
            data_smooth_max_roughness.obs_SAR_plane = obs_SAR_plane_max_roughness;
            data_smooth_max_roughness.Q_SAR = Q_SAR;
        else
            temp_greens = load(model_greens_path);
            data_smooth_max_roughness.Gss_SAR = temp_greens.Gss_SAR;
            data_smooth_max_roughness.Gds_SAR = temp_greens.Gds_SAR;
            data_smooth_max_roughness.ll_SAR = temp_greens.ll_SAR;
        end



         % plotting_model(data_smooth_max_roughness,figprop,'Roughest')

    end



    %% Search for the case with maximum smoothness distribution
    if number_files-skip_files>=1  && smoothness_exclude~=1  % Smoothness was used as variable

        % the number of the run within the selected file
        ix_run = ix_total_max_smoothness-(ix_file_max_smoothness-skip_files-1)*n_runs_file;

        % loading the datafile
        outdata_temp = load([outdata_path_files '_' num2str(ix_file_max_smoothness) '.mat']);
        m_opt_smooth = outdata_temp.m_keep(:,ix_run);

        % design matrix including greens coefficients and plane
        if GPS_alone~=1
            A_plane = outdata_temp.A_design(:,end-2:end);
            n_sar = length(obs_SAR);
            A_plane_SAR = A_plane(1:n_sar,:);

            if sigma_exclude~=1
                xslope_plane=m_opt_smooth(end-4);
                yslope_plane=m_opt_smooth(end-3);
                offset_plane=m_opt_smooth(end-2);
            else
                xslope_plane=m_opt_smooth(end-3);
                yslope_plane=m_opt_smooth(end-2);
                offset_plane=m_opt_smooth(end-1);
            end

            % Computing the fitted plane to save as output
            obs_SAR_plane_smooth = A_plane_SAR*[xslope_plane ;yslope_plane; offset_plane];
        end

        % convert rake slip to dip and strike slip and then ivert using greens
        % to the slip.
        rake_smooth=m_opt_smooth(1:n_patch);
        rake_slip_smooth = m_opt_smooth(n_patch+1:2*n_patch);
        strike_slip_smooth =m_opt_smooth(n_patch+1:2*n_patch).*cos(rake_smooth*pi/180);
        dip_slip_smooth = m_opt_smooth(n_patch+1:2*n_patch).*sin(rake_smooth*pi/180);
        m_slip_smooth(1:n_patch)=strike_slip_smooth;                                    % strike slip
        m_slip_smooth(n_patch+1:2*n_patch)=dip_slip_smooth;                             % dip slip
        if GPS_alone~=1
            m_slip_smooth(2*n_patch+1:2*n_patch+3)=m_opt_smooth(2*n_patch+1:2*n_patch+3);   % plane
        end
        clear rake_smooth 


        U = outdata_temp.A_design*(m_slip_smooth)';
        % splitting surface deformations in a GPS and InSAR component
        if GPS_alone==1
            obs_GPS_opt_smooth = U;
        else
            obs_SAR_opt_smooth = U(1:n_sar);    
            obs_GPS_opt_smooth = U(n_sar+1:end);
        end

        clear U

        % Decomposing GPS modelled observations vector into multiple components 
        obs_GPS_opt_smooth_E = obs_GPS_opt_smooth(1:n_GPS);
        obs_GPS_opt_smooth_N = obs_GPS_opt_smooth(n_GPS+1:2*n_GPS);
        obs_GPS_opt_smooth_U = obs_GPS_opt_smooth(2*n_GPS+1:end);

        % give the strike and dip greens parameters of the InSAR such the model
        % can be shown.
        data_smooth.dip_slip = dip_slip_smooth;
        data_smooth.strike_slip = strike_slip_smooth;
        data_smooth.rake_slip= rake_slip_smooth;
        data_smooth.m =  outdata_temp.model_interface;
        data_smooth.ll_GPS = ll_GPS;
        data_smooth.cal0 = cal0;  
        data_smooth.transitions = transitions;
        data_smooth.GPS_alone=GPS_alone;
        data_smooth.obs_GPS_E = obs_GPS_E;
        data_smooth.obs_GPS_N = obs_GPS_N;
        data_smooth.obs_GPS_U = obs_GPS_U;
        data_smooth.obs_GPS_opt_E = obs_GPS_opt_smooth_E;
        data_smooth.obs_GPS_opt_N = obs_GPS_opt_smooth_N;
        data_smooth.obs_GPS_opt_U = obs_GPS_opt_smooth_U;
        data_smooth.A_design =  outdata_temp.A_design;
        data_smooth.Q_GPS = Q_GPS;
        data_smooth.outdata_path_files = outdata_path_files;
        data_smooth.GPS_LOS_struct = outdata.obs_GPS_LOS;
        data_smooth.checkerboard = checkerboard;

        if GPS_alone~=1
            data_smooth.ll_SAR = ll_SAR;
            data_smooth.obs_SAR = obs_SAR;
            data_smooth.obs_SAR_opt = obs_SAR_opt_smooth;
            data_smooth.obs_SAR_plane = obs_SAR_plane_smooth;
            data_smooth.Q_SAR = Q_SAR;
        else
            temp_greens = load(model_greens_path);
            data_smooth.Gss_SAR = temp_greens.Gss_SAR;
            data_smooth.Gds_SAR = temp_greens.Gds_SAR;
            data_smooth.ll_SAR = temp_greens.ll_SAR;
        end



         % plotting_model(data_smooth,figprop,'MaxSmoothness')

    end

end


%%




patch_n = 50;
% plotting for patch n





if GPS_alone==1
    
    
    
    
    
    h1 = figure('name',[' dislocation patch ' num2str(ix_big_patches)],'position',[5         502        2552         420]);
    subplot(3,1,1)
    plot(rake_slip_keep_all(ix_big_patches,:),'k-')
    title('rake slip')
    subplot(3,1,2)
    plot(rake_keep_all(ix_big_patches,:),'k-')
    title('rake')
    subplot(3,1,3)
    plot(hyper_keep,'k-')
    title('smoothness')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'disloc_' num2str(ix_big_patches) '.eps'])
    end


    h1 = figure('name',[' dislocation patch ' num2str(patch_n)],'position',[5          42        2552         420]);
    subplot(3,1,1)
    plot(rake_slip_keep_all(patch_n,:),'k-')
    title('rake slip')
    subplot(3,1,2)
    plot(rake_keep_all(patch_n,:),'k-')
    title('rake')
    subplot(3,1,3)
    plot(hyper_keep,'k-')
    title('smoothness')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'disloc_' num2str(patch_n) '.eps'])
    end
    
    
    h1 = figure('name',[' P Post'],'position',[5          42        2552         420]);
    plot(P_post_keep_all,'k-')
    title('P Post')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'P_keep.eps'])
    end
        
    

    
else
    h1 = figure('name',[ 'dislocation patch ' num2str(ix_big_patches)],'position',[5         502        2552         420]);
    subplot(3,1,1)
    plot(rake_slip_keep_all(ix_big_patches,:),'k-')
    title('rake slip')
    subplot(3,1,2)
    plot(rake_keep_all(ix_big_patches,:),'k-')
    title('rake')
    subplot(3,1,3)
    plot(hyper_keep,'k-')
    title('smoothness')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'disloc_' num2str(ix_big_patches) 'a.eps'])
    end

    h1 = figure('name',['dislocation patch ' num2str(ix_big_patches)],'position',[5         502        2552         420]);
    subplot(3,1,1)
    plot(plane_keep(1,:),'k-')
    title('slope X')
    subplot(3,1,2)
    plot(plane_keep(2,:),'k-')
    title('slope Y')
    subplot(3,1,3)
    plot(plane_keep(3,:),'k-')
    title('offset')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'disloc_' num2str(ix_big_patches) 'b.eps'])
    end
    
    h1 = figure('name',[ 'dislocation patch ' num2str(patch_n)],'position',[5          42        2552         420]);
    subplot(3,1,1)
    plot(rake_slip_keep_all(patch_n,:),'k-')
    title('rake slip')
    subplot(3,1,2)
    plot(rake_keep_all(patch_n,:),'k-')
    title('rake')
    subplot(3,1,3)
    plot(hyper_keep,'k-')
    title('smoothness')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'disloc_' num2str(patch_n) 'a.eps'])
    end
    
    h1 = figure('name',['dislocation patch ' num2str(patch_n)],'position',[5          42        2552         420]);
    subplot(3,1,1)    
    plot(plane_keep(1,:),'k-')
    title('slope X')
    subplot(3,1,2)
    plot(plane_keep(2,:),'k-')
    title('slope Y')
    subplot(3,1,3)
    plot(plane_keep(3,:),'k-')
    title('offset')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'disloc_' num2str(patch_n) 'b.eps'])
    end
    
    
        
    h1 = figure('name',[' P Post'],'position',[5          42        2552         420]);
    plot(P_post_keep_all,'k-')
    title('P Post')
    set(gcf,'PaperPositionMode','auto')
    if save_fig_flag ==1
        print(h1,'-depsc','-r150',[outdata_path_files filesep 'P_keep.eps'])
    end
        
    
end











