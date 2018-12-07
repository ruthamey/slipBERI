function [x,y]=circharm(coeffs,phi,plot_flag)
%CIRCHARM circular harmonics
%   [x,y]=circharm(coeffs,phi,plot_flag)
%   coeffs = strength of each harmonic
%   phi = rotation of each harmonic (excluding first)
%
%   Andy Hooper Dec 2010

if nargin<3
   plot_flag=1;
end

theta=linspace(0,2*pi,101);
theta=theta(1:100);
%r=(1+coeffs(1)/sqrt(2*pi))*ones(1,100);
%coeffs=coeffs(2:end);

P=legendre_poly(length(coeffs)-1,cos(theta));
for i=2:size(P,1)
    P(i,:)=circshift(P(i,:),[0,round(phi(i-1)/2/pi*100)]);
end
r=sum(P.*repmat(coeffs(:),1,100));
%plot(theta,P)
%figure

%for i=2:2:length(coeffs)
%     r=r+coeffs(i)*cos(i/2*theta)/sqrt(pi);
%end
%
%for i=3:2:length(coeffs)
%     r=r+coeffs(i)*sin((i-1)/2*theta)/sqrt(pi);
%end

x=r.*cos(theta);
y=r.*sin(theta);

if plot_flag~=0
    n=ceil(sqrt(length(coeffs)));
    subplot(n,n,1)
    plot(x,y,'b.')
    title('\Sigma C_lY_l(\theta)')
    axis equal
    for i=2:length(coeffs)
        subplot(n,n,i)
        x1=P(i,:).^2*coeffs(i).*cos(theta);
        y1=P(i,:).^2*coeffs(i).*sin(theta);
        plot(x1(P(i,:)>=0),y1(P(i,:)>=0),'b.')
        hold on
        plot(x1(P(i,:)<0),y1(P(i,:)<0),'r.')
        hold off
        title(['l=',num2str(i)])
        axis equal
    end
end
