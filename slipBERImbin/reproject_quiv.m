function [xquivmag,yquivmag,zquivmag] = reproject_quiv(quiv_mags,strcol,dipcol)
% REPROJECT_QUIV: Reproject quiver arrows from along-strike / down-dip to
% X,Y,Z coordinate system.
% from Tom Ingleby, May 2017

xquivmag = quiv_mags(:,1).*sind(strcol) - quiv_mags(:,2).*cosd(dipcol).*cosd(strcol);
yquivmag = quiv_mags(:,1).*cosd(strcol) + quiv_mags(:,2).*cosd(dipcol).*sind(strcol);
zquivmag = quiv_mags(:,2).*sind(dipcol);

end

