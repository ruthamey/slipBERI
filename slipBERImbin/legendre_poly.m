function [P]=legendre_poly(n,x)
%LEGENDRE_POLY return Legendre polyonmoials
%   function [P]=legendre_poly(n,x)
%  
%
%   Andy Hooper Dec 2010


P=ones(n+1,length(x));
P(2,:)=x;

for i=1:n-1
     P(i+2,:)=((2*i+1)*x.*P(i+1,:)-i*P(i,:))/(i+1);
end
    
