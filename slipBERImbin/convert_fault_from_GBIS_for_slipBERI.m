function [] = convert_fault_from_GBIS_for_slipBERI(inversion_result_file, n_along_strike_patches, n_down_dip_patches, outname)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Load results from GBIS
load(inversion_result_file, 'invResults')
optimalmodel = cell2mat(invResults.optimalmodel(1));
faultlength = optimalmodel(1);
faultwidth = optimalmodel(2);
faultdepth = optimalmodel(3);
dip = abs(optimalmodel(4));
strike = optimalmodel(5)-180;   % GBIS uses a different strike/dip convention to slipBERI
xcenter = optimalmodel(6);
ycenter = optimalmodel(7);
strikeslip = optimalmodel(8);
dipslip = optimalmodel(9);

% Calculate rake
rake = tand(dipslip/strikeslip);

% Convert x and y from local coordinate system to lat long
load(inversion_result_file, 'geo')
origin = geo.referencePoint;
llh=local2llh([xcenter; ycenter]/1000,origin);
xcenter_ll = llh(1);
ycenter_ll = llh(2);

% Calculate bottom depth
bottomdepth = sind(abs(dip))*faultwidth;

     %strike %dip %rake  %x center   % y center  % length        %topdepth %bottom depth, n_along_strike_patches, n_down_dip_patches, smoothing, 1
M = [strike, dip, rake, xcenter_ll, ycenter_ll, faultlength/1000,  0,     bottomdepth/1000,   n_along_strike_patches, n_down_dip_patches      1,     1];
dlmwrite(outname,M, 'delimiter', ' ')

end

