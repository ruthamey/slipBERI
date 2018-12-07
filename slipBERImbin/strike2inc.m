function [xinc,yinc]=strike2inc(strike)

% function [xinc,yinc]=strike2inc(strike)
% 
% Calculates the x and y increments (assuming a magnitude of 1)
% of a vector with a given strike direction.
%
% A bit brute force. Agricultural, if you will.
% 
% gjf, 16-jun-2004
% 
% genius in action (tm)

deg2rad=pi/180;

if ((strike>=0) & (strike<90))

	str=strike*deg2rad;
	xinc=sin(str);
	yinc=cos(str);

else if ((strike>=90) & (strike<180))

	str=(strike-90)*deg2rad;
	xinc=cos(str);
	yinc=-1*sin(str);

else if ((strike>=180) & (strike<270))

	str=(strike-180)*deg2rad;
	xinc=-1*sin(str);
	yinc=-1*cos(str);

else if ((strike>=270) & (strike<360))

	str=(strike-270)*deg2rad;
	xinc=-1*cos(str);
	yinc=sin(str);

end
end
end 
end
