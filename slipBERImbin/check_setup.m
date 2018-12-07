% A quick code that plots the data and fault, to check that everything is
% in the correct place and there have been no units mix-ups
%
% R.M.J.Amey 2018


figure;

for i= 1: max(InSAR_identifyer)

    subplot(1,max(InSAR_identifyer),i)
    hold on;
    if strcmp(data.InSAR_datafile, 'none') ~= 1
        scatter( locs_InSAR(1,InSAR_identifyer==i), locs_InSAR(2,InSAR_identifyer==i), 100, d_InSAR(InSAR_identifyer==i), 'filled');   
    end

    if strcmp(data.atolls_datafile, 'none') ~= 1
        scatter( locs_atolls(1,:), locs_atolls(2,:), 100, d_atolls);
    end

    if strcmp(data.GPS_datafile_2d, 'none') == 0
        quiver( locs_GPS_2d(1,1:2:end), locs_GPS_2d(2,1:2:end), d_GPS_e_2d, d_GPS_n_2d, 'k');
    end
    
    if strcmp(data.GPS_datafile_3d, 'none') == 0
        quiver( locs_GPS_3d(1,1:3:end), locs_GPS_3d(2,1:3:end), d_GPS_e_3d, d_GPS_n_3d, 'k');
    end

        for j = 1 : n_fault_strands
                x = [fault_coords(j,1), fault_coords(j,3)];
                y = [fault_coords(j,2), fault_coords(j,4)];
                hold on;
                plot(x, y, 'm', 'Linewidth', 2)
        end

        axis equal
        colormap('jet')

        xlabel('UTM x')
        ylabel('UTM y')
        ylabel(colorbar, 'LOS displacement (m)')

end

if strcmp(testing.testing_mode, 'yes') == 1
    
    load(testing.making_model, 'u');
    from_making_model = load(testing.making_model, 'locs');
    scale = 340; 
    
    % true data
    displ_actual_E = u(1,:);
    displ_actual_N = u(2,:);
    quiver( from_making_model.locs(1,:), from_making_model.locs(2,:), displ_actual_E* scale, displ_actual_N* scale, 'b', 'Autoscale', 'off');
    hold on;
    
        for j = 1 : n_fault_strands
                x = [fault_coords(j,1), fault_coords(j,3)];
                y = [fault_coords(j,2), fault_coords(j,4)];
                hold on;
                plot(x, y, 'm', 'Linewidth', 2)
        end

        axis equal
        colormap('jet')

        xlabel('UTM x')
        ylabel('UTM y')
    
end