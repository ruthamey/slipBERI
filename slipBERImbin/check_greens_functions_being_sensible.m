% A hastingly created, quite possibly very specific code that plots the
% forward model of the current slip values, to see if the mechanism is
% broadly correct
%
% R.M.J.Amey 2018

    slip_curr = m_curr(m_identifyer==1);
    
%     if strcmp(invert.solve_for_InSAR_ramp, 'yes') == 1 || strcmp(invert.solve_for_InSAR_offset, 'yes') == 1
%         offset_curr = offset_or_ramp_temp;
%     else
%         offset_curr = zeros(length(d),1);
%     end
    
    figure('position', [100, 350, 1800, 600])
    ax1 = subplot(1,3,1);
    title('Data');
    
    % First, plot InSAR, as long as there is some
    if strcmp(data.InSAR_datafile, 'none') ~= 1
        
        scatter( locs_InSAR(1,:), locs_InSAR(2,:), [], d_InSAR, 'filled')
  
        
        axis equal
        colorbar
        %xlabel('UTM x, I think')
        %ylabel('UTM y')
    end
    
    scale_factor = 50000;
    
    % Then plot GPS, as long as there is some
    if strcmp(data.GPS_datafile_2d, 'none') == 0 || strcmp(data.GPS_datafile_3d, 'none') == 0
           subplot(1,3,1); hold on;
           quiver( locs_GPS_unique(1,:), locs_GPS_unique(2,:), d_GPS_e*scale_factor, d_GPS_n*scale_factor, 'k', 'Autoscale', 'off');   % i.e. plot every third row of the locations against E and N...
           ylabel('UTM y')
           title('Data');    
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
        


    % Second, plot surface displacement predicted from slip solution
    ax2 = subplot(1,3,2);
    
%    if strcmp(invert.variable_rake, 'yes') == 1
        G = G_curr;
%     else
%         G = G;
%     end  
    
    d_hat =  (G * slip_curr) + offset_or_ramp_temp;
    
    subplot(1,3,2); hold on;
    title('Data');
    
    if strcmp(data.InSAR_datafile, 'none') ~= 1    
        scatter( locs_InSAR(1,:), locs_InSAR(2,:), [], d_hat(1:length(locs_InSAR)), 'filled')
        hold on;
    end
    
    if strcmp(data.GPS_datafile_2d, 'none') + strcmp(data.GPS_datafile_3d, 'none') == 0
        d_hat_e_2d = d_hat((length(locs_InSAR)+1):2:(length(locs_InSAR) + length(d_GPS_2d)));
        d_hat_n_2d = d_hat((length(locs_InSAR)+2):2:(length(locs_InSAR) +length(d_GPS_2d)));
        d_hat_e_3d = d_hat(((length(locs_InSAR) + length(locs_GPS_2d)+1)):3:end);
        d_hat_n_3d = d_hat(((length(locs_InSAR) + length(locs_GPS_2d)+2)):3:end);
        d_hat_e = [d_hat_e_2d; d_hat_e_3d];
        d_hat_n = [d_hat_n_2d; d_hat_n_3d];
    elseif strcmp(data.GPS_datafile_3d, 'none') == 0;   % just 3d 
        d_hat_e = d_hat((length(locs_InSAR)+1):3:end);
        d_hat_n = d_hat((length(locs_InSAR)+2):3:end);
    elseif strcmp(data.GPS_datafile_2d, 'none') == 0;   % just 2d 
        d_hat_e = d_hat((length(locs_InSAR)+1):2:end);
        d_hat_n = d_hat((length(locs_InSAR)+2):2:end);
    end
    
    if strcmp(data.GPS_datafile, 'yes') == 1;
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

    title('Model (MAP)');

    
    % Third, plot residuals
    ax3 = subplot(1,3,3);
    residuals = (d - d_hat);%.^2;        % the first half are InSAR residuals, if we have InSAR. The second half are GPS residuals, if we have GPS.
    
    if strcmp(data.InSAR_datafile, 'none') ~= 1     % if we have InSAR data, then the first (length(locs_InSAR)) entries of d are for InSAR
        scatter( locs_InSAR(1,:), locs_InSAR(2,:), [], residuals(1:length(locs_InSAR)), 'filled');
        hold on;
    end
    
    if strcmp(data.atolls_datafile, 'none') ~= 1     % if we have InSAR data, then the first (length(locs_InSAR)) entries of d are for InSAR
        scatter( locs_atolls(1,:), locs_atolls(2,:), [], residuals((length(d_InSAR)+length(d_GPS)+1):end), 'filled');
        hold on;
    end
  
    
    if strcmp(data.GPS_datafile_2d, 'none') + strcmp(data.GPS_datafile_3d, 'none') == 0
           subplot(1,3,3); hold on;
           % 2d
           quiver( locs_GPS(1, 1:2: length(locs_GPS_2d)), locs_GPS(2, 1:2:length(locs_GPS_2d)), residuals((length(locs_InSAR)+1):2:(length(locs_InSAR) + length(d_GPS_2d)))'*scale_factor, residuals((length(locs_InSAR)+2):2:(length(locs_InSAR) +length(d_GPS_2d)))'*scale_factor, 'k', 'Autoscale', 'off');   % 2D so plot every second row of the locations against E and N...
           % then 3d
           quiver( locs_GPS(1, (length(locs_GPS_2d)+1):3:end), locs_GPS(2, (length(locs_GPS_2d)+2):3:end), d_hat(((length(locs_InSAR) + length(locs_GPS_2d)+1)):3:end)'*scale_factor, d_hat(((length(locs_InSAR) + length(locs_GPS_2d)+2)):3:end)'*scale_factor, 'k', 'Autoscale', 'off');   % i.e. 3D so plot every third row of the locations against E and N...        
           title('Data'); 
    elseif strcmp(data.GPS_datafile_2d, 'none') == 0; 
           subplot(1,3,3); hold on;
           quiver( locs_GPS(1,1:2: length(d_GPS)), locs_GPS(2,1:2: length(d_GPS)), residuals((length(locs_InSAR)+1):2:end)'*scale_factor, residuals((length(locs_InSAR)+2):2:end)'*scale_factor, 'k', 'Autoscale', 'off');   % i.e. plot every third row of the locations against E and N...
           title('Data'); 
    elseif strcmp(data.GPS_datafile_3d, 'none') == 0;
           subplot(1,3,3); hold on;
           quiver( locs_GPS(1,1:3: length(d_GPS)), locs_GPS(2,1:3: length(d_GPS)), residuals((length(locs_InSAR)+1):3:end)'*scale_factor, residuals((length(locs_InSAR)+2):3:end)'*scale_factor, 'k', 'Autoscale', 'off');   % i.e. plot every third row of the locations against E and N...
           title('Data');      
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

    title('Residuals');

    % Make sure all colorbars are same scale, only display on RHS
    
    c = max(abs([min(d_InSAR), max(d_InSAR)])); % Calculate maximu value for symmetric colormap
    caxis([-c c])
    
    
    %MaxC = max( [max(d_InSAR), max(d_hat(1:length(locs_InSAR))), max(residuals(1:length(locs_InSAR)))]);  % max colour is max InSAR value - ignore GPS arrows
    %MinC = min( [min(d_InSAR), min(d_hat(1:length(locs_InSAR))), min(residuals(1:length(locs_InSAR)))]);  % max colour is max InSAR value - ignore GPS arrows
    %subplot(1,3,1); caxis([ MinC MaxC]); colorbar('off');
    %subplot(1,3,2); caxis([ MinC MaxC]); colorbar('off');
    %subplot(1,3,3); caxis([ MinC MaxC]);
    subplot(1,3,1); caxis([ -c c]); colorbar('off');
    subplot(1,3,2); caxis([ -c c]); colorbar('off');
    subplot(1,3,3); caxis([ -c c]);
    if strcmp(data.GPS_datafile, 'none') ~= 1
        ylabel(colorbar, 'LOS displacement (m)')    % note that GPS has been projected into LOS
    end
    %colormap('hsv')
    [redbluecmap] = redblue;
    colormap(flipud(redbluecmap))
    if strcmp(data.atolls_datafile, 'none') ~= 1;
        ylabel(colorbar, 'vertical displacement (m)')    % note that GPS has been projected into LOS
    end
    
    linkaxes([ax1, ax2, ax3], 'xy')