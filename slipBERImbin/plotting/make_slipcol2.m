
for r = 1:401
col(r,1) = 1;
  if (r<134)
     col(r,2) = 1;
     col(r,3) = 0 + (133-(r-1))/133;
  elseif (r<267)
     col(r,3)=0 ;
     col(r,2) = 0 + (133-(r-134))/133;  
  elseif (r<=401)
     col(r,2)=0 ;
     col(r,3) = 1 - (136-(r-267))/136;
  end
end
