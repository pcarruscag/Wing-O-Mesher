%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [X,Y,M,r] = Inflation(Xin,Yin,data)

N = data.N;  Min = data.M;
l1 = data.firstLayer;  r = data.rate;

printf('  Generating inflation layer... '); fflush(stdout);
tic

% determine line length along const qsi lines
etaL = zeros(Min,N+1);
for i=2:Min
  etaL(i,:) = etaL(i-1,:)+sqrt((Xin(i,:)-Xin(i-1,:)).^2+...
                               (Yin(i,:)-Yin(i-1,:)).^2);
  
endfor

% use the shortest line to determine the required new number of points in eta
% and their distribution (results are more reasonable this way)
[L,idx] = min(etaL(Min,:));
etaR = sqrt(Xin(:,idx).^2+Yin(:,idx).^2);

si = [0 l1];
while (L-si(end))/(si(end)-si(end-1)) > 0.5
  M = numel(si);
  % constant rate
  ds = r*(si(M)-si(M-1));
  % limit aspect ratio (based on perimeter of ellipse with a/b = 10)
  dt = 1.3*pi*interp1(etaL(:,idx),etaR,si(end),'spline','extrap')/(N-1);
  ds = min(ds,dt*data.nearAR);
  % append to si
  si = [si si(end)+ds];
endwhile
si(end) = L;
si(end-1) = 0.5*(si(end)+si(end-2));
si /= L;
% new number of eta points
M++;

% re-mesh along constant qsi lines with forced orthogonality
X = zeros(M,N+1); Y = zeros(M,N+1);
m = data.numOrthog;
if m < 0
  w = zeros(1,M);
else
  w = exp(log(0.01)*[zeros(1,m) 1:M-m]/(M-m)); w(end) = 0;
endif

for i=2:N
  s = etaL(:,i);  x = Xin(:,i);  y = Yin(:,i);
  
  if i<N
    theta = atan2(Xin(1,i+1)-Xin(1,i-1), Yin(1,i-1)-Yin(1,i+1));
  else
    theta = atan2(Xin(1,3)-Xin(1,1), Yin(1,1)-Yin(1,3));
  endif
  
  xn = x(1)+si*L*cos(theta);
  yn = y(1)+si*L*sin(theta);
  
  leDist = sqrt((x(1)+0.5*data.chord)^2 + y(1)^2)/(0.005*data.t*data.chord);
  teDist = sqrt((x(1)-0.5*data.chord)^2 + y(1)^2)/(0.005*data.t*data.chord);
  off = 1 - max(1-teDist,0).^4.*(1+4*teDist) - max(1-leDist,0).^4.*(1+4*leDist);
  
  X(:,i) = (1-off*w).*interp1(s,x,si*s(end),'spline','extrap') + off*w.*xn;
  Y(:,i) = (1-off*w).*interp1(s,y,si*s(end),'spline','extrap') + off*w.*yn;
endfor

X(:, 1 ) = X(:,N); Y(:, 1 ) = Y(:,N);
X(:,N+1) = X(:,2); Y(:,N+1) = Y(:,2);

printf('%.2fs\n',toc); fflush(stdout);

endfunction
