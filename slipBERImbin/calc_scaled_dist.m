function [ scaled_dists, anglebetween ] = calc_scaled_dist( n_fault_strands_for_smoothing, disloc_model, a_as, a_dd, first_patch_in_strand_for_smoothing, last_patch_in_strand_for_smoothing, along_strike_sep_dist, n_along_strike_patches, n_down_dip_patches, fault_strand_togetherness)

% calc_scaled_dist gives a matrix of svaled distances between two slip patches
% scaled by the von karman scaling relations ax and az (Mai & Beroza 2002)
%
% It returns a symmetric n x n matrix giving the distance between slip
% patches, denoted by the row number. e.g. row 2, column 3 gives the
% distance between slip patch 2 and slip patch 3.
%
% e.g. if there were a 2 x 3 matrix like this:
%    (1 3 5
%     2 4 6)
%
% Where the along-strike separation distance is 'a', and the down-dip
% separation distance is 'b'
%
% This function would return the following symmetric matrix:
%
% (         0
% (         a                          0
% (         b                  sqrt(a^2 + b^2)          0
% (    sqrt(a^2 + b^2)                 a                b                  0
% (         2a               sqrt(2a^2 + b^2)           a            sqrt(a^2 + b^2)    0 
% (  sqrt(2a^2 + b^2)              2a            sqrt(a^2 + b^2)           a            b       0 )     
%
%
% Where the diagonals are equal to 0, because the distance between a slip
% patch and itself is 0
%
% It also returns the angle between the centre of slip patches IN RADIANS, 
% to use later when calculating the direction dependent Hurst number
%
% Ruth Amey 18-dec-2014
% rmja      12-feb-2015  Solving for correlation lengths, multiple fault strands
% rmja      3-jan-2016   Adapted from calc_scaled_dist so that it can cope with changing number of slip patches along strike

% Housekeeping
scaled_x_dists = [];
scaled_z_dists = [];
scaled_x_dists_temp = [];
scaled_z_dists_temp = [];
scaled_dists = [];
anglebetween = [];
fault_im_on = 1;

% Work out each fault as if it's a separate fault
for i = 1: n_fault_strands_for_smoothing

        x_locs = [];
            n_fault_strands_for_smoothing_on_each_strand = sum(fault_strand_togetherness==i);
                starter_value = 0;
                for n = fault_im_on:(fault_im_on+n_fault_strands_for_smoothing_on_each_strand-1)
                    x_locs_temp = along_strike_sep_dist(n)*(1: (n_along_strike_patches(n))) - 1/2*along_strike_sep_dist(n) + starter_value;
                    starter_value = x_locs_temp(end) + (1/2*along_strike_sep_dist(n));
                    x_locs_temp = repmat(x_locs_temp,n_down_dip_patches(n), 1);
                    x_locs_temp = reshape(x_locs_temp, [],1);
                    x_locs = [x_locs; x_locs_temp];
                    x_locs_temp = [];
                end
             fault_im_on = fault_im_on + n_fault_strands_for_smoothing_on_each_strand;

    z_locs_temp = (disloc_model(9, first_patch_in_strand_for_smoothing(i):last_patch_in_strand_for_smoothing(i)) + disloc_model(8,first_patch_in_strand_for_smoothing(i):last_patch_in_strand_for_smoothing(i))) / 2;  % using depth from disloc model (have to calculate centre from top depth and bottom depth)

    % Calculate the distances between the centres of each fault patch on the same fault strand
    x_dists_temp =  pdist2(x_locs,x_locs);
    z_dists_temp =  pdist2(z_locs_temp',z_locs_temp');
    
    % Scale by von karman ax, az
    scaled_x_dists_temp(:,:) = x_dists_temp / a_as(i);
    scaled_z_dists_temp(:,:) = z_dists_temp / a_dd(i);
    
    % Work out total scaled distance between each patch on the same fault strand and angle between them
    scaled_dist_temp(:,:) = sqrt( scaled_x_dists_temp(:,:).^2 + scaled_z_dists_temp(:,:).^2);
    angle_temp(:,:) = atan(scaled_z_dists_temp(:,:)./scaled_x_dists_temp(:,:));
   
    % Make one big matrix of distances - distances between patches on different fault strands will be 0
     scaled_dists = blkdiag(scaled_dists, scaled_dist_temp); % since we only need to know the distances between patches on the SAME fault strand, for vk sigma_s calculation
     anglebetween = blkdiag(anglebetween, angle_temp);
     scaled_x_dists = blkdiag(scaled_x_dists, scaled_x_dists_temp);
     scaled_z_dists = blkdiag(scaled_z_dists, scaled_z_dists_temp);
     
     scaled_x_dists_temp = [];
     scaled_z_dists_temp = [];
     x_dists_temp = [];
     z_dists_temp = [];
     angle_temp = [];
     scaled_dist_temp = [];
end
    
% Set diagonals to 0 rather than NaN
anglebetween(isnan(angle_temp)) = 0;

end

