function [ux,uy,uz,err] = dc3d3(alpha,x,y,dip,al1,al2,aw1,aw2,disl1,disl2) 
%% function [ux,uy,uz,err] = dc3d3(alpha,x,y,dip,al1,al2,aw1,aw2,disl1,disl2)
%%
%% Direct Translation of Pete Clarke's dc3d3.c into octave code
%% dc3d3.c is itself a translation of Okada's dc3d.f subroutine into C code
%% ...which has been pruned to compute ux,uy,uz for z=0 only, no tensile 
%% components, depth=0
%%
%% tjw 11-sept-02

%% Constants
     F0 = 0;
     F1 = 1;
     F2 = 2;
     PI2 = 2*pi;
     EPS = 1.0e-6;

     [woo,nx]=size(x);

%% Other vector gubbins

     disl1_3=[disl1; disl1; disl1];
     disl2_3=[disl2; disl2; disl2];

     u = zeros(3,nx);
     dub = zeros(3,nx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                % 
%%dccon0 "subroutine"
% Calculates medium and fault dip constants
     c0_alp3 = (F1-alpha)/alpha;

     pl8 = PI2/360;
     c0_sd = sin(dip*pl8);
     c0_cd = cos(dip*pl8);
     if (abs(c0_cd) < EPS)
       c0_cd = F0;
       if (c0_sd>F0)  c0_sd = F1; end
       if (c0_sd<F0)  c0_sd = -F1; end
     end
     c0_cdcd = c0_cd*c0_cd;
     c0_sdcd = c0_sd*c0_cd;
%                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     p=y*c0_cd;
     q=y*c0_sd;

     jxi = (((x-al1).*(x-al2))<=F0);
     jet = (((p-aw1).*(p-aw2))<=F0);
     
     for k = 1:2
       if (k==1) et = p-aw1; else et = p-aw2; end 
       for j = 1:2
         if (j==1) xi = x-al1; else xi = x-al2; end


% 

%       disp([x(1) y(1) p(1) q(1) jxi(1) jet(1)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                %
%%dccon2 "subroutine"
% calculates station geometry constants for finite source

       etq_max=max([abs(et); abs(q)]); 
       dc_max = max([abs(xi); etq_max]);
 
%       dc_max = max(abs(xi),max(abs(et),abs(q)));
%       if ((abs(xi/dc_max) < EPS) || (abs(xi) < EPS)) xi = F0; end
%       if ((abs(et/dc_max) < EPS) || (abs(et) < EPS)) et = F0; end
%       if ((abs(q/dc_max) < EPS) || (abs(q) < EPS)) q = F0; end


       xi=xi-(abs(xi./dc_max) < EPS).*xi;
       xi=xi-(abs(xi) < EPS).*xi;
       et=et-(abs(et./dc_max) < EPS).*et;
       et=et-(abs(et) < EPS).*et;
       q=q-(abs(q./dc_max) < EPS).*q;
       q=q-(abs(q) < EPS).*q;


       dc_xi = xi; dc_et = et; dc_q = q;


 %      disp([dc_max(1) xi(1) et(1) q(1)])


       c2_r = sqrt(dc_xi.^2 + dc_et.^2 + dc_q.^2);

       if (sum((c2_r) == F0) > 0)  
           disp(['singularity error']);
           ux=0;
           uy=0;
           uz=0;
           err=1; 
           return;
       end

       c2_y = dc_et*c0_cd + dc_q*c0_sd;
       c2_d = dc_et*c0_sd - dc_q*c0_cd;

%       if (dc_q == F0) c2_tt = F0; else c2_tt = atan(dc_xi.*dc_et./(dc_q.*c2_r)); end

       c2_tt = atan(dc_xi.*dc_et./(dc_q.*c2_r));
       c2_tt = c2_tt-c2_tt.*(dc_q == F0);

%       if ((dc_xi < F0) && (dc_q == F0) && (dc_et == F0)) c2_x11 = F0; 
%       else 
%         rxi = c2_r + dc_xi;
%         c2_x11 = F1/(c2_r*rxi);
%       end

       c2_x11 = (ones(1,nx)*F1)./(c2_r.*(c2_r+dc_xi));
       c2_x11 = c2_x11-((dc_xi < F0).*(dc_q == 0).*(dc_et == 0)).*c2_x11;

%       if ((dc_et < F0) && (dc_q == F0) && (dc_xi ==F0))
%         if ((c2_r-dc_et) < 1e-14) disp(['dccon2 a']);disp([c2_3-dc_et, c2_r, dc_et, dc_q, dc_xi]); end
%         c2_ale = -log(c2_r-dc_et);
%         c2_y11 = F0 ;
%       else
%	 ret = c2_r + dc_et;
%         if (ret < 1e-14) disp(['dccon2 b']); disp([ret,c2_r,dc_et,dc_q,dc_xi]); end
%         c2_ale = log(ret);
%         c2_y11 = F1/(c2_r*ret);
%       end


       c2_ale_msk1 = ((dc_et < F0) .* (dc_q == F0) .* (dc_xi ==F0));
       c2_ale_msk2 = 1 - c2_ale_msk1;

       ret1 = c2_r-dc_et;
       ret2 = c2_r+dc_et;

       c2_ale = (c2_ale_msk1.*-log(ret1)) + (c2_ale_msk2.*log(ret2));
       c2_y11 = (ones(size(ret2))*F1)./(c2_r.*ret2);

%       disp([c2_r(1) c2_y(1) c2_d(1) c2_tt(1) c2_x11(1) c2_y11(1) c2_ale(1)])

                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%         if ( ( (q==F0) && (((jxi==1) && (et==F0)) || ((jet==1)&&(xi==F0))) ) || (c2_r==F0) ) 
%           ux = 0; uy = 0; uz = 0;
%	   err = 2; % singular problems: return error code
%	   return;
%         end

         if ( sum( (q == F0) .* ( ((jxi == 1).*(et == F0)) + ((jet == 1).*(xi == F0)) ) + (c2_r == F0) ) > 0)
           ux = 0; 
           uy = 0;
           uz = 0;
	   err = 2; % singular problems: return error code
	   return;
         end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                %
%% ub "subroutine"
%  part B of displacement and strain at depth due to buried fauls in semi-infinite medium

       rd = c2_r + c2_d;

       if ( sum(rd < 1e-14) > 0) 
         disp(['ub']);disp([rd ,c2_r, c2_d,xi,et,q]);end


       if (c0_cd ~= F0)

%	 if (xi==F0) ai4 = F0; 
%         else
%           xx=sqrt(xi.^2+q.^2); % xx replaces x in original subroutine
%           ai4 = F1/c0_cdcd * (xi/rd*c0_sdcd + F2*atan((et.*(xx+q*c0_cd) + xx*(c2_r+xx)*c0_sd) / (xi*(c2_r+xx)*c0_cd)) );
%         end
%         ai3 = (c2_y*c0_cd/rd - c2_ale + c0_sd*log(rd)) / c0_cdcd;

        xx=sqrt(xi.^2+q.^2);

        ai4 = (xi ~= 0).* (  (F1/c0_cdcd) * (xi./rd*c0_sdcd + F2*atan((et.*(xx+q*c0_cd) +...
             xx.*(c2_r+xx)*c0_sd) ./ (xi.*(c2_r+xx)*c0_cd)) )  );
        ai3 = (c2_y*c0_cd./rd - c2_ale + c0_sd*log(rd)) / c0_cdcd;

       else

	 rd2 = rd.*rd;
	 ai3 = (et./rd + c2_y.*q./rd2 - c2_ale) / F2; %%% HERE
         ai4 = (xi.*c2_y./rd2)/F2;

       end

       ai1 = -xi./rd*c0_cd - ai4*c0_sd;
       ai2 = log(rd) + ai3*c0_sd;

       qx = q.*c2_x11;
       qy = q.*c2_y11;

%       disp([ai1(1) ai2(1) ai3(1) ai4(1)])


%strike-slip contribution
       if (disl1 ~=0) 
         du2(1,:) = - xi.*qy  - c2_tt  - c0_alp3*ai1*c0_sd;
	 du2(2,:) = - q./c2_r          + c0_alp3*c2_y./rd*c0_sd;
	 du2(3,:) =   q.*qy            - c0_alp3*ai2*c0_sd;

%         size(disl1)
%         size(du2)
 
	 dub = (disl1_3.*du2)/PI2; 
       else
	 dub = zeros(3,nx);
       end
       
%       disp([du2(1,1) du2(2,1) du2(3,1)]) 


%dip-slip contribution
       if (disl2 ~=F0)
 	 du2(1,:) = - q./c2_r          + c0_alp3*ai3*c0_sdcd;
	 du2(2,:) = - et.*qx - c2_tt   - c0_alp3*xi./rd*c0_sdcd;
	 du2(3,:) =   q.*qx            + c0_alp3*ai4*c0_sdcd;
	 dub = dub + (disl2_3/PI2).*du2;
       end

%       disp([du2(1,1) du2(2,1) du2(3,1)]) 


%                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	 du(1,:) = dub(1,:);
         du(2,:) = dub(2,:)*c0_cd - dub(3,:)*c0_sd;
         du(3,:) = dub(2,:)*c0_sd + dub(3,:)*c0_cd;
         if((j+k)~=3) u = u+du; else u = u-du;end
       end
     end

%       disp([dub(1,1) dub(2,1) dub(3,1)]) 



     ux = u(1,:); 
     uy = u(2,:); 
     uz = u(3,:); 
     err = 0;

%       disp([ux(1) uy(1) uz(1)]) 
