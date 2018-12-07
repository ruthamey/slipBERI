function [  ] = doplot3d( faults, colorbar )

%%% Script to plot fault geometry in 3D
%
% Borrowed from slipinv (Funning 2005)

[rp,cp] = size(faults);
d2r=pi/180;
hold off

for p = 1:rp

xcent = faults(p,1);
ycent = faults(p,2);
length = faults(p,7);
strike = faults(p,3)*d2r;
top = faults(p,8);
bot = faults(p,9);
slip = faults(p,6);
dip = faults(p,4)*d2r;
%dip = 45*d2r;

x_surf1 = xcent - 0.5*length*sin(strike);
x_surf2 = xcent + 0.5*length*sin(strike);
y_surf1 = ycent - 0.5*length*cos(strike);
y_surf2 = ycent + 0.5*length*cos(strike);

x_top1 = x_surf1 + top*cos(strike)/tan(dip);
x_top2 = x_surf2 + top*cos(strike)/tan(dip);
y_top1 = y_surf1 - top*sin(strike)/tan(dip);
y_top2 = y_surf2 - top*sin(strike)/tan(dip);


x_bot1 = x_surf1 + bot*cos(strike)/tan(dip);
x_bot2 = x_surf2 + bot*cos(strike)/tan(dip);
y_bot1 = y_surf1 - bot*sin(strike)/tan(dip);
y_bot2 = y_surf2 - bot*sin(strike)/tan(dip);

X(1:4,p) = [x_top1 x_top2 x_bot2 x_bot1]'/1000;
Y(1:4,p) = [y_top1 y_top2 y_bot2 y_bot1]'/1000;
Z(1:4,p) = [-top  -top  -bot  -bot]'/1000;


fill3(X(:,p),Y(:,p),Z(:,p),slip)
axis equal
make_slipcol2
colormap(col);
hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

caxis([0 max(faults(:,6))])
colormap(colorbar);

end
