function [disloc_model, spatial_model1, spatial_model2, spatial_model3, spatial_model4, seg_coords, oksar_out, alongstrike, downdip, fault_length_meters, fault_width_meters, n_fault_strands, strike, fault_segment_togetherness, twicesmoothed_togetherness] =  process_faultdata_centre(fault, invert, use_local_coordinate_system, origin, testing)

% function [disloc_model, spatial_model1, spatial_model2, spatial_model3, fault_coords, faultseg_coords, oksar_out] = 
%   process_faultdata (input_filename)
%
% This is a function which formats fault model data, from the input
% file specified, into a form that the Okada function 'disloc3d'
% can read. It also generates a spatial model which relates the 
% various fault patches spatially, and a dataset that will be used 
% to generate 'oksar'-friendly output later on.
%
% outputs (according to rmja):
% disloc_model - one column per patch with details in
%       disloc_model(1) = x centre of patch;
%       disloc_model(2)= y centre of patch;
%       disloc_model(3)= strike of patch;
%       disloc_model(4)= dip of patch;
%       disloc_model(5)= rake of patch;
%       disloc_model(6)= slip of patch, set initially to 1;
%       disloc_model(7)= patch_length
%       disloc_model(8,l)=top(i)+(k-1)*vertwidth(i) = top of patch;
%       disloc_model(9,l)=disloc_model(8,l)+vertwidth(i) = bottom of patch;
% spatial_model1 - slip patch number (labelled down each column)
% spatial_model2 - slip patch width
% spatial_model3 - slip patch height
% spatial_model4 - identifyer. I don't use this, but I think this is if you want to do layered earth model, so different elastic parameters with depth
% seg_coords - coordinates of fault. fault segment end coordinates x,y,x,y. One row per fault strand.
%
% nicked from gfj, hacked by rmja

% input data file
%  input_filename
  fin = fopen(fault.fault_descriptor_file);
  fault_input = fscanf(fin,'%lf',[12,inf]);
  fault_input = fault_input';
  fclose(fin);
  
  
% RUTH HACKS **************************************************************
n_fault_strands = size(fault_input, 1);

if n_fault_strands > 1
    
   % when I plot them flat using imagesc, want to plot LHS of fault being NORTH or WEST. so depending on strike, might need to do some flipping, depending on the order at which you look at the fault strands 
   
   % if the fault is striking approximately southwards...
   % we want to look at the most northerly fault first, if fault is striking > 90 and < 180. I.e. if it's striking approx south, we want to plot FIRST patch in the TOP LEFT, so this is fine.
   % so we need to make sure the faults are given to us in the order: most northerly --> least northerly  
%    if any(fault_input(:,1) > 90 & any(fault_input(:,1)< 180));
%        
%        [~, order ] = sort( fault_input(:,5), 'descend'); 
%        if all(diff(order)<0) == 1; % if it doesn't start with the most northerly one (we want the order to be [1 2 3], so the diff(order) is always > 0). if not we need to do some switchin'
%         fault_input = fault_input(order,:);
%         disp('have switched fault strand order in processs_faultdata_ceter, so most northerly fault is first, so that imagesc plots properly') 
%        end
%        
%    % but if the fault is striking approximately northwards...
%    % want to consider the top LEFT of all the faults first, and then flip them all at the end. so make sure we consider the most southern fault first.
%     
%    elseif any(fault_input(:,1) > 180 | any(fault_input(:,1) < 90));
%    
%        [~, order ] = sort( fault_input(:,5), 'descend'); 
%        if all(diff(order)>0) == 1; % if it doesn't start with the most northerly one (we want the order to be [1 2 3], so the diff(order) is always > 0). if not we need to do some switchin'
%         fault_input = fault_input( flipud(order),:);
%         disp('have switched fault strand order in processs_faultdata_ceter, so most southerly fault is first, will sort out later, I promise') 
%        end
%    
%    
%    end
  
   
   
end






  
% set up variables

strike=fault_input(:,1) ;
dip=fault_input(:,2) ;
rake=fault_input(:,3) ;


% RUTH HACKS **************************************************************

xc=fault_input(:,4) ;
yc=fault_input(:,5) ;

if strcmp(testing.testing_mode, 'no') == 1;
    if xc < 0
        xc = xc+360;
    end

    if strcmp(use_local_coordinate_system, 'yes') == 1;

        llh = [xc'; yc'; zeros(1,length(yc))];
         xy = llh2local(llh, origin);
         xc = (xy(1,:)*1000)';
         yc = (xy(2,:)*1000)'; 

    end

    if strcmp(fault.fault_coordinate_unit, 'long/lat') == 1;
         % CONVERT TO UTMX and UTMY if not given in that format
         %disp('check in process_faultdata_centre_ruthhack that you''re doing it right in situations when you DON''T convert to local coordinate system.')
         [xc, yc, ~ ] = ll2utm(yc, xc);
    end
    
else   % if in testing mode
     xc=fault_input(:,4)*1000 ;
     yc=fault_input(:,5)*1000 ; 
    %disp('RUTH you commented out * 1000 in process_faultdata_cetner')
end

fault_length=fault_input(:,6) ;          % ruth-hack - change 'length' to 'fault_length'
fault_length_meters = fault_length * 1000 ;      % ruth-hack - CONVERT TO METERS
%**************************************************************************
    


top=fault_input(:,7) ;
bottom=fault_input(:,8) ;
fault_width = abs(top - bottom);         % ruth-hack - I want to know the width of the fault, too
fault_width_meters = fault_width * 1000;        % ruth-hack - CONVERT TO METERS
alongstrike=fault_input(:,9) ;
downdip=fault_input(:,10) ;
%affil=fault_input(:,11);      % ruth-hack commented out. i'm usingaffiliation in a different way ( i want one value per fault strand, notone value per patch )
affil = 1;
fault_segment_togetherness = fault_input(:,11);
twicesmoothed_togetherness = fault_input(:,12);


%seg_coords = fault_input(:,1:4);

deg2rad=pi/180 ;

% calculate number of model parameters

[m,woo] = size(xc);  %number of fault 'segments'

% work out dimensions of the spatial models, and set 'em up
%rws=max(downdip);                   % ruth hack changed from 'rws=max(downdip);' since sometimes we have faults that go below the surface so we need to sum the downdip patches, and not sum the alongstrikes
%clms=sum(alongstrike);
if any(top~=0)    % if any of the faults start below ground then we need to sum all the down dip patches, and they'll all have the same number of along strike patches
    clms = max(alongstrike);       % this is necessary because sometimes we have multiple strands laterally (e.g. napa) and sometimes we have them going down depth (e.g. nepal)
    rws = sum(downdip);
else    % if all the faults start at the surface but there are a few extending laterally then we need to sum all the alongstrike patches, and they'll all have the same number of downdip patches
    clms = sum(alongstrike);
    rws = downdip(1);
end
spatial_model1=zeros(rws,clms);
spatial_model2=zeros(rws,clms);
spatial_model3=zeros(rws,clms);
spatial_model4=ones(rws,clms)*-1;


% set count variables to zero

l=0 ;
p=0 ;

% set reversal flag to zero

reversed=0 ;

% set up loops in which to calculate the desired parameters
% (e.g. strike, dip, length, width...) and put them into the
% rows of a matrix 'disloc_model' which is in the correct
% format for the Okada routine 'disloc3d'


% loop through entries in the input data

for i=1:m ;

% check the dip of the fault plane... if it is greater than 90
% degrees, take appropriate action - switching end point
% coordinates, adjusting the dip, and reversing the direction of
% strike-slip (this is to allow reversals in dip along strike)

    if (dip(i)>90)
 
% reversal flag        
        
        reversed=1 ;
        
% reverse strike

        strike(i)=strike(i)+180 ;
        if (strike(i)>=360)
		strike(i)=strike(i)-360;
	end
        
% dip adjustment:

       dip(i)=180-dip(i) ;   
       
% rake reversal:        
        
       rake(i)=-1*rake(i) ;
    
    elseif (dip(i)<=90)

        
% reversal flag        
        
        reversed=0 ;
        
    end
    
% calculate fault patch lengths

    kmlength(i)=fault_length(i)/alongstrike(i) ;
    patch_length(i)=kmlength(i)*1000 ;                   % ruth hack, rename to patch_length
    
    
% calculate fault patch widths

    width(i)=(faultwidth(dip(i), top(i), bottom(i)))/downdip(i); 
    %vertlength = (bottom(i)-top(i))*1000 ;
    %width(i) = (vertlength/sin(dip*(pi/180)))/downdip(i) ;
    
    vertwidth(i) = (bottom(i)-top(i))*1000/downdip(i);

    
   
% calculate strike-slip and dip-slip components
    
    rrake(i)=(rake(i)+90)*deg2rad ;
    dipslip(i)=cos(rrake(i)) ;
    strslip(i)=sin(rrake(i)) ;
    
% calculate the x and y increments

    [xinctmp,yinctmp]=strike2inc(strike(i));
    xinc(i)=xinctmp*patch_length(i);
    yinc(i)=yinctmp*patch_length(i);

% calculate the fault segment end coordinates
    
    seg_coords(i,1)=xc(i)-(xinctmp*patch_length(i)*alongstrike(i)/2);
    seg_coords(i,2)=yc(i)-(yinctmp*patch_length(i)*alongstrike(i)/2);
    seg_coords(i,3)=xc(i)+(xinctmp*patch_length(i)*alongstrike(i)/2);   
    seg_coords(i,4)=yc(i)+(yinctmp*patch_length(i)*alongstrike(i)/2);


% loop through along-strike subdivisions   


% if reversed, loop through coordinates backwards

    if (reversed==1)
     startj = alongstrike(i)-(alongstrike(i)+1)/2;
     incj   = -1;
     endj   = 1-(alongstrike(i)+1)/2;
    elseif (reversed==0)
     startj = 1-(alongstrike(i)+1)/2;
     incj   = 1;
     endj   = alongstrike(i)-(alongstrike(i)+1)/2;
    end


    for j=startj:incj:endj 

% increment along-strike division counter

        p=p+1 ;
        
% calculate end coordinates for surface trace of the subdivided
% fault patches:       
        
%        tempx1(i,j)=((x2(i)-x1(i))*(j-1)/alongstrike(i))+x1(i) ;
%        tempy1(i,j)=((y2(i)-y1(i))*(j-1)/alongstrike(i))+y1(i) ;
%        tempx2(i,j)=((x2(i)-x1(i))*j/alongstrike(i))+x1(i) ;
%        tempy2(i,j)=((y2(i)-y1(i))*j/alongstrike(i))+y1(i) ;
        
%       tempxcentre(i,j)=tempx1(i,j)+(tempx2(i,j)-tempx1(i,j))/2 ;
%        tempycentre(i,j)=tempy1(i,j)+(tempy2(i,j)-tempy1(i,j))/2 ;

        tempxcentre=xc(i)+j*xinc(i);
        tempycentre=yc(i)+j*yinc(i) ;
        
% dip-number counter            
            
if any(top~=0)   % if any of the fault is buried
     dip_row_counter_start = [0; (cumsum(downdip))];
     dip_row_counter_start(end) =[];
     dip_column_counter = reshape( repmat(1:alongstrike(1), sum(downdip), 1)', numel(repmat(1:alongstrike(1), sum(downdip), 1)), 1);
else
     dip_row_counter_start = zeros( m, 1);
     dip_column_counter = 1:10000;         % abritrarily large...
end        
        
              
% loop through down-dip subdivisions        
        
        for k=1:downdip(i)
            
% increment fault-element counter            
            
            l=l+1 ;       
            
            
% dip counter           
            
           dip_row_counter = dip_row_counter_start(i);
            
% disloc_model output:            
 

% x location (centre of projection of fault plane)
            disloc_model(1,l)=tempxcentre;

% y location (centre of projection of fault plane)
            disloc_model(2,l)=tempycentre;

% strike
            disloc_model(3,l)=strike(i);
           
% dip 
            disloc_model(4,l)=dip(i);

% rake
            disloc_model(5,l)=rake(i);

% slip
            disloc_model(6,l)=1;            

% length
            disloc_model(7,l)=patch_length(i);      % ruth-hack - change to patch length

% hmin      
            disloc_model(8,l)=top(i)*1000+(k-1)*vertwidth(i);      % ruth-hack, used to be    disloc_model(8,l)=top(i)+(k-1)*vertwidth(i);
% hmax
            disloc_model(9,l)=disloc_model(8,l)+vertwidth(i);

% spatial_model output: % NB octave only deals with 2D arrays not 3D           
            
% set fault element number (layer 1)

            spatial_model1(dip_row_counter + k,dip_column_counter(p))=l ;
            
% set fault element length (layer 2)            
            
            spatial_model2(dip_row_counter + k,dip_column_counter(p))=patch_length(i) ;      % ruth-hack - change to patch length
            
% set fault element width (layer 3)     
            
            spatial_model3(dip_row_counter + k,dip_column_counter(p))=width(i);

% set smoothing 'affiliation' (layer 4)

            %spatial_model4(dip_row_counter + k,dip_column_counter(p))=affil(i);        % ruth-hack          
            
% oksar_out output:

% set x coordinate (row 1)

            oksar_out(1,l)=tempxcentre ;

% set y coordinate (row 2)            
            
            oksar_out(2,l)=tempycentre;
         
% set strike (row 3)            
            
            oksar_out(3,l)=strike(i) ;
            
% set dip (row 4)
            
            oksar_out(4,l)=dip(i) ;
            
% set rake (row 5) 
            
            oksar_out(5,l)=rake(i) ;
            
% set length (row 6)

            oksar_out(6,l)=kmlength(i) ;
            
% set top (row 7)
            
            oksar_out(7,l)=top(i)+((k-1)*...
                (bottom(i)-top(i))/downdip(i)) ;
 
% set bottom (row 8)
            
            oksar_out(8,l)=bottom(i)-((downdip(i)-k)*...
                (bottom(i)-top(i))/downdip(i)) ;
            
            
        end
        

    end

end
    
% RUTH HACKS **************************************************************
n_fault_strands = size(fault_input, 1);

if any(top~=0)
    
    % this totally won't work so don't bother
    
else

    if n_fault_strands > 1

        
        disp('have commented out the end of process_faultdata_centre_ruthhack so PLOTTING FLAT WILL BE BACKWARDS but hopefully the rest will work better');

%        if any(fault_input(:,1) > 180 | any(fault_input(:,1) < 90)); 
% 
%            disp('now flipping disloc_model, so that when we plot in imagesc, number ONE patch in every fault strand is the top LEFT. so the LEFT of our flat fault is NORTH.')
% 
%            order = flipud(spatial_model1);
%            %order = fliplr(order);
%            %order = reshape(order, 1, sum(alongstrike)*bottom(1));
%            order = reshape(order, 1, sum(alongstrike)*downdip(1));
%            order = order';
%            disloc_model = disloc_model(:, order);
%            disloc_model = fliplr(disloc_model);
%            %spatial_model1 = flip(ud(fliplr(disloc_model));     % THIS IS FINE - this is just numbering the patches, and so we want to keep it as it is .
%            spatial_model2 = flipud(fliplr(spatial_model2));     % THIS IS NOT QUITE - this assumes all patches on the same fault strand have the same property
%            spatial_model3 = flipud(fliplr(spatial_model3));     % THIS IS NOT QUITE - this assumes all patches on the same fault strand have the same property
%            spatial_model4 = flipud(fliplr(spatial_model4));     % THIS IS NOT QUITE - this assumes all patches on the same fault strand have the same property
% 
%        end
    end

end

%**************************************************************************




end

    
