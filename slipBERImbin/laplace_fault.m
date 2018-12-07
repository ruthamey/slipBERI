function [L]=laplace_fault(n,m,r,free_edge)
%LAPLACE_FAULT return laplacian design matrix for n x m matrix
%   L=laplace(n,m,r,free_edge)
%   n is number of rows or cols in first dimension
%   m is number of rows or cols in second dimension
%   r is ratio of patch length in 1st dim over length in 2nd dim
%   free_edge is 0 if fault is buried, 
%                1 if edge at beginning of first dimension is free
%                2 if edge at end of first dimension is free 
%
%   Example (numbers give order of patches):
%
%           m=4
%       1  4  7  10  
%   n=3 2  5  8  11
%       3  6  9  12
%
%   free_edge=1 means top row is free, 2 means bottom row is free. 
%
%   Andy Hooper, 2011
%
%==============================================================================
% 09/2012 AH Corners of free edge set correctly
% 10/2014 AH renamed from laplace.m to avoid conflict with matlab function
%==============================================================================

if nargin<3
   r=1; % ratio of patch length in 1st dim over length in 2nd dim
end

if nargin<4
   free_edge=0; % free edge, 
end

if free_edge>2
    error('free_edge must be 0, 1 or 2')
end

r2=r^2;
c=-(2*r2+2); % value in middle (4 in case of square patches)

%left col
L1=zeros(n-2,n*m);
L1(1:n-2,1:n-2)=eye(n-2);
L1(1:n-2,2:n-1)=L1(1:n-2,2:n-1)+eye(n-2)*c;
L1(1:n-2,3:n)=L1(1:n-2,3:n)+eye(n-2);
L141=L1;
L1(1:n-2,n+2:2*n-1)=L1(1:n-2,n+2:2*n-1)+eye(n-2)*r2;
L=L1;

%main body
Lm=zeros(n,n*m);
Lm(1:n,1:n)=eye(n)*r2;
Lm(2:n,n+1:2*n-1)=eye(n-1);
Lm(1:n,n+1:2*n)=Lm(1:n,n+1:2*n)+eye(n)*c;
if free_edge==2
    Lm(end,2*n)=c+r2;
end
if free_edge==1
    Lm(1,n+1)=c+r2;
end
Lm(1:n-1,n+2:2*n)=Lm(1:n-1,n+2:2*n)+eye(n-1);
Lm(1:n,2*n+1:3*n)=eye(n)*r2;
L=[L;Lm];

for i=1:m-3
   L=[L;circshift(Lm,[0,i*n])];
end 

L2=circshift(L141,[0,(m-1)*n]);
L2(1:n-2,(m-2)*n+2:(m-1)*n-1)=L2(1:n-2,(m-2)*n+2:(m-1)*n-1)+eye(n-2)*r2;

L=[L;L2];

%corners
LC=zeros(4,n*m);

if free_edge==1
    LC(1,1:2)=[c+r2,1];
    LC(3,(m-1)*n+1:(m-1)*n+2)=[c+r2,1];
else
    LC(1,1:2)=[c,1];
    LC(3,(m-1)*n+1:(m-1)*n+2)=[c,1];
end
if free_edge==2
    LC(2,n-1:n)=[1,c+r2];
    LC(4,end-1:end)=[1,c+r2];
else
    LC(2,n-1:n)=[1,c];
    LC(4,end-1:end)=[1,c];
end
LC(1,n+1)=r2;
LC(2,2*n)=r2;
LC(3,(m-2)*n+1)=r2;
LC(4,end-n)=r2;

L=[L;LC];
