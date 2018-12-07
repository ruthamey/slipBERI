function [ x,y,u,v ] = prepare_plot_rake( rake, slip, spatial_model2, spatial_model3, n_down_dip_patches, n_along_strike_patches, disloc_model )
% prepare_plot_rake takes the output of slipBERI and transforms the value
% of rake and slip into an output that can be plotted using the quiver plot
%
% This function takes the value of rake and slip for each patch, and
% changes this polar-coordinate style measurement into (u,v), which is the
% magnitude of the x and y values of your arrow
%
% It then reads the centre of each patch from disloc model, and leaves you
% with something nice to plot
%
% rmja  14-aug-2015  rainy friday afternoon
%
%
% turn rake into u and v
%theta = 180 - rake;               % turn rake into a polar coordinate angle
% theta = rake-180;
% theta = degtorad(theta);
% [u,v] = pol2cart(theta, slip);           % slip gives the magnitude of the arrow
%
%
% find centre of each patch, so you can plot it there
% x = disloc_model(1,:);
% y = disloc_model(2,:);
%
% SO, at the moment I'm just plotting using imagesc, which means that the
% patches on row one, column one have centre point (1,1). so I just need to
% work out which row all the patches are on and that can be my centre
% points
% if strcmp(testing.are_you_testing, 'no') == 1;
%     theta = 180-rake;
%     theta = degtorad(theta);
%     [u,v] = pol2cart(theta, slip');           % slip gives the magnitude of the arrow
%   x = disloc_model(1,:);
%   y = disloc_model(2,:);
%   
% elseif strcmp(testing.are_you_tesiing, 'yes') == 1;

slip_rows = size(slip,1);
slip_columns = size(slip,2);

if slip_columns > slip_rows % if rows are longer than columns
   slip = slip';             % make it a column vector 
end

    theta = rake-180;
    theta = deg2rad(theta);
        if size(theta,1) ~= size(slip,1)  % transpose if necessary, so they're both column vectors
        theta = theta';
        end
    [u,v] = pol2cart(theta, slip);           % slip gives the magnitude of the arrow
    
%     N = length(slip);
%     slip_patches = (1:N).';
%     fault_shape = reshape( slip_patches, n_down_dip_patches, n_along_strike_patches);
%     for i = 1 : N
%        [row, column] = find(fault_shape == i); % since the patches are numberered by order    
%        y(i, 1) = row;               % column vector
%        x(i, 1) = column;            % column vector
%     end
    
     row = repmat([(1:n_down_dip_patches)'], [1, n_along_strike_patches]);
    %column = repmat([(1:n_along_strike_patches)'], [1,n_down_dip_patches])';
    
%     y = spatial_model2 .* row ;          % dd sep distance * row number
%     y = reshape( y, n_down_dip_patches*n_along_strike_patches,1)./1000;
     %y = ((disloc_model(9,:) - ((disloc_model(9,:) -disloc_model(8,:)) / 2))/ 1000)';        % down-dip centre of each fault patch
    %y = spatial_model3(1)/1000/sin(deg2rad(disloc_model(4,1)))/2 + ( spatial_model3(1)/1000/sin(deg2rad(disloc_model(4,1))) * (row-1))
      y = spatial_model3(1)/2/1000 + ((row-1) * spatial_model3(1)/1000);
    y = reshape( y, n_down_dip_patches*n_along_strike_patches,1);
     
%     x = spatial_model3 .* column- (spatial_model3/2);       % as sep distance * column number    
      x = (cumsum(spatial_model2,2) - (spatial_model2/2))/ 1000;
      x = reshape( x, n_down_dip_patches*n_along_strike_patches,1);
   
%end


end

