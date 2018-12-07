function [U,flag]=disloc3d2(m,x,lambda,mu)
% DISLOC3D    [U,flag]=disloc3d(m,x,lambda,mu)
%
% Returns the deformation at point 'x', given dislocation
% model 'm'.  lambda and mu are the Lame parameters
%
% Both 'm' and 'x' can be matrices, holding different models
% and observation coordinates in the columns.  In this case,
% the function returns the deformation at the points specified
% in the columns of 'x' from the sum of all the models in the
% columns of 'm'.  'x' must be 2xi (i = number of observation
% coordinates) and 'm' must be 9xj (j = number of models).
%
% The coordinate system is as follows: east = positive X,
% north = positive Y. Calculations assume Z=0.
%
% The main output is 'U', the three displacement components:
% east, north, and up (on the rows)).  This output has
% the same number of columns as 'x'.
%
% Output 'flag' is set for a singularity.
%
% The dislocation model is specified as: east,north,strike,dip,
% rake,slip,length,hmin,hmax. The coordinates (east,north)
% specify the midpoint of the surface projection of the fault
% plane.
%
% Modified to use cut down, octave version of dc3d.f: dc3d3.m.
% No Strains, No Tensile Components, Depth of observations = 0.
%
% Vectorised for optimal speed in MATLAB 17-may-2005

     DEG2RAD = (2*pi/360);
     nfaults = size(m,2);
     nx      = size(x,2);

%initiate matrix for output displacements, U
     U=zeros(3,nx);

% calculate alpha
     alpha = (lambda+mu)/(lambda+2*mu);

%loop over models
     for i=1:nfaults
%convert model into correct inputs for dc3d3.m
       flt_x = ones(1,nx)*m(1,i);
       flt_y = ones(1,nx)*m(2,i);
       strike = m(3,i);
       dip = m(4,i);
       rake = m(5,i);
       slip = m(6,i);
       length = m(7,i);
       hmin = m(8,i);
       hmax = m(9,i);

       rrake = (rake+90)*DEG2RAD;
       sindip = sin(dip*DEG2RAD);
       w = (hmax-hmin)/sindip;
       ud = ones(1,nx)*slip*cos(rrake);
       us = ones(1,nx)*-slip*sin(rrake);
       halflen = length/2;
       al2 = ones(1,nx)*halflen;
       al1 = -al2;
       aw1 = ones(1,nx)*hmin/sindip;
       aw2 = ones(1,nx)*hmax/sindip;
%reject data which breaks the surface
       if (sum(hmin < 0)>0) disp(['ERROR: Fault top above ground surface']);end;
       hmin=hmin+(hmin == 0)*0.00001;

       sstrike = (strike+90)*DEG2RAD;
       ct = cos(sstrike); 
       st = sin(sstrike);

%loop over points

         X=ct*(-flt_x+x(1,:))-st*(-flt_y+x(2,:));
         Y=ct*(-flt_y+x(2,:))+st*(-flt_x+x(1,:));


         [ux,uy,uz,err]=dc3d4(alpha,X,Y,-dip,al1,al2,aw1,aw2,us,ud);


         if (err~=0) disp(['error code in dc3d3... stopping!']);flag = err,return;else flag=0;end
	 U(1,:) = U(1,:) + ct*ux + st*uy;
         U(2,:) = U(2,:) -st*ux + ct*uy;
         U(3,:) = U(3,:) + uz;

     end
