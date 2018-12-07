function width = faultwidth(dip, top, bottom)

% function width = faultwidth(dip, top, bottom)
%
% faultwidth calculates the down-dip width of a fault (in metres)
% given its dip (in decimal degrees), and the depths to 
% the top and bottom of the fault (in kilometres)
%

vertlength = (bottom-top)*1000 ;

width = vertlength/sin(dip*(pi/180)) ;
