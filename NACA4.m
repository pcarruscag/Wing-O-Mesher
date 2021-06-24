%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [X,Y] = NACA4(m,p,t,N,k,w)
%NACA4 surface coordinates of a 4 digit NACA airfoil
% m - camber
% p - location of maximum camber
% t - thickness
% N - number of points
% k - factor to bias seeding towards TE
% w - factor to control the minimum TE size

u = 0:1/(N/2-1):1;
u = (1-k)*u + k*(1-(u-1).^2);
x = (1-cos(w*pi*u))*0.5;
x /= x(end);

% camber line
yc=zeros(size(x));
yc(x<0.1*p)=m/p^2*x(x<0.1*p).*(0.2*p-x(x<0.1*p));
yc(x>=0.1*p)=0.01*m*(1-x(x>=0.1*p))/(1-0.1*p)^2.*(1+x(x>=0.1*p)-0.2*p);

% thickness distribution
yt=0.05*t*(0.2969*sqrt(x)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);

% camber line derivative
dyc=zeros(size(x));
dyc(x<0.1*p)=2*m/p^2*(0.1*p-x(x<0.1*p));
dyc(x>=0.1*p)=2*m/(10-p)^2*(0.1*p-x(x>=0.1*p));

% suction and pressure side points
xu=x-yt.*sin(atan(dyc));
xl=x+yt.*sin(atan(dyc));
yu=yc+yt.*cos(atan(dyc));
yl=yc-yt.*cos(atan(dyc));

% round the TE
tte = 2*0.05*t*0.0021;
keep = x<(1-tte); % these points we do not modify
idx = [1:sum(keep) numel(x)+1 N+2-fliplr(1:sum(keep))];
idx_i = [1:numel(x) N+2-fliplr(1:numel(x))];
X=[xu(keep) 1.0 fliplr(xl(keep))];
Y=[yu(keep) 0.0 fliplr(yl(keep))];

% final geometry
X = interp1(idx,X,idx_i,'spline')-0.5;
Y = interp1(idx,Y,idx_i,'spline');

end
