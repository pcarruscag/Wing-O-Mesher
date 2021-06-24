%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
nLE = data.nTh1;
nTAR = data.nTh2;

function [Xq0,Yq0,Xq1,Yq1] = semiCircleMesh(refangle,cutratio,x,y,le=false)
  % mesh a semi-circle shape using 2 blocks of quadrilaterals

  % compute angles to decide where to cut
  theta = atan2(diff(y),diff(x));
  for idx=1:numel(theta)
    if theta(idx) < (pi/2-refangle)
      theta = pi/2-theta(idx);
      idx++;
      break
    endif
  endfor

  cut1 = idx;

  % the length of the LE segment
  lcut = 2*sum(sqrt(diff(x(1:idx)).^2+diff(y(1:idx)).^2));

  % based on it we decide where to cut next
  l = [0 cumsum(sqrt(diff(x).^2+diff(y).^2))];
  while l(idx) < cutratio*lcut
    idx++;
  endwhile

  cut2 = idx;

  % locations of key points
  % first cut
  xu1 = x(cut1);
  yu1 = y(cut1);

  % second cut
  xu2 = x(cut2);
  yu2 = y(cut2);

  % bottom points
  cut3 = numel(x)-cut1+1;
  cut4 = numel(x)-cut2+1;

  % first cut
  xl1 = x(cut3);
  yl1 = y(cut3);

  % second cut
  xl2 = x(cut4);
  yl2 = y(cut4);

  % triple point (bifurcation)
  xu3 = xu1+yu1*cos(theta)/(1+sin(theta));
  yu3 = yu1/(1+sin(theta));

  % the lower side angle here may be different from the reference
  theta = atan((y(cut3+1)-y(cut3))/(x(cut3+1)-x(cut3)))+pi/2;
  xl3 = xl1-yl1*cos(theta)/(1+sin(theta));
  yl3 = yl1/(1+sin(theta));

  % end of second block
  % draw a line from xu2 to xl2, cut such that ||P1-P3|| = ||P2-P4||
  w = 1 - sqrt((xu3-xu1)^2+(yu3-yu1)^2) / sqrt((xu2-xl2)^2+(yu2-yl2)^2);
  xu4 = xu2*w + xl2*(1-w);
  yu4 = yu2*w + yl2*(1-w);

  xl4 = xl2*w + xu2*(1-w);
  yl4 = yl2*w + yu2*(1-w);

  % first we mesh the inner block as that will determine the size
  % of the cells in the other 3
  nH = 2*(cut1-1)+1;
  nW = cut2-cut1+1;
  dW = sqrt((xu4-xu3)^2+(yu4-yu3)^2)/nW;
  if le
    nC = max(floor((0.5*(xu3+xl3)-x(1))/dW),3);
  else
    nC = max(ceil((nH/(1-2*(1-w))-nH)/2),3);
  endif
  nH2 = nH+2*nC-2;

  Xq1 = zeros(nH,nW);  Yq1 = zeros(nH,nW);

  Xq1(1,:)   = interp1([1 nW], [xl3 xl4], 1:nW);
  Xq1(end,:) = interp1([1 nW], [xu3 xu4], 1:nW);
  Yq1(1,:)   = interp1([1 nW], [yl3 yl4], 1:nW);
  Yq1(end,:) = interp1([1 nW], [yu3 yu4], 1:nW);

  w = (0:nH-1)'/(nH-1);
  Xq1 = w*Xq1(end,:) + (1-w)*Xq1(1,:);
  Yq1 = w*Yq1(end,:) + (1-w)*Yq1(1,:);

  % upper block
  Xu = zeros(nC,nW);  Yu = zeros(nC,nW);
  Xu(1,:)   = Xq1(end,:);
  Xu(end,:) = x(cut1:cut2);
  Yu(1,:)   = Yq1(end,:);
  Yu(end,:) = y(cut1:cut2);

  w = (0:nC-1)'/(nC-1);
  Xu = w*Xu(end,:) + (1-w)*Xu(1,:);
  Yu = w*Yu(end,:) + (1-w)*Yu(1,:);

  % lower block
  Xl = zeros(nC,nW);  Yl = zeros(nC,nW);
  Xl(1,:)   = fliplr(x(cut4:cut3));
  Xl(end,:) = Xq1(1,:);
  Yl(1,:)   = fliplr(y(cut4:cut3));
  Yl(end,:) = Yq1(1,:);

  w = (0:nC-1)'/(nC-1);
  Xl = w*Xl(end,:) + (1-w)*Xl(1,:);
  Yl = w*Yl(end,:) + (1-w)*Yl(1,:);

  % tip block
  Xq0 = zeros(nH,nC);  Yq0 = zeros(nH,nC);
  Xq0(:,1)   = [x(cut3:end)'; x(2:cut1)'];
  Xq0(:,end) = Xq1(:,1);
  Yq0(:,1)   = [y(cut3:end)'; y(2:cut1)'];
  Yq0(:,end) = Yq1(:,1);

  w = (0:nC-1)/(nC-1);
  Xq0 = Xq0(:,end)*w + Xq0(:,1)*(1-w);
  Yq0 = Yq0(:,end)*w + Yq0(:,1)*(1-w);

  % merge blocks
  Xq1 = [Xl; Xq1(2:end-1,:); Xu];
  Yq1 = [Yl; Yq1(2:end-1,:); Yu];
  
  % make last section uniform
  [M,N] = size(Xq1);
  Xq1(:,end) = interp1([1 M], [Xq1(1,end) Xq1(end,end)], 1:M, "linear")';
  Yq1(:,end) = interp1([1 M], [Yq1(1,end) Yq1(end,end)], 1:M, "linear")';

endfunction


### Leading and Trailing Edge Semi-Circles ###

x = X(1,2:end);
y = Y(1,2:end);

[Xq0,Yq0,Xq1,Yq1] = semiCircleMesh(data.refAngle,1.5,x,y,true);

cut1 = ceil(rows(Xq0)/2)+columns(Xq1)-1;
cut2 = data.N-cut1+1;

x2 = -[0.5*(x(data.N/2)+x(data.N/2+1))...
       fliplr(x(2:data.N/2))...
       fliplr(x(data.N/2+1:end-1))];
x2(end+1) = x2(1);
y2  = [0.5*(y(data.N/2)+y(data.N/2+1))...
       fliplr(y(2:data.N/2))...
       fliplr(y(data.N/2+1:end-1))];
y2(end+1) = y2(1);

[Xq5,Yq5,Xq4,Yq4] = semiCircleMesh(data.refAngle,1.0,x2,y2);
clear x2 y2

% remove construction point and flip back
[mid,~] = size(Xq4);  mid = ceil(mid/2);
Xq4 = -fliplr([Xq4(1:mid-1,:);  Xq4(mid+1:end,:)]);
Yq4 =  fliplr([Yq4(1:mid-1,:);  Yq4(mid+1:end,:)]);

[mid,~] = size(Xq5);  mid = ceil(mid/2);
Xq5 = -fliplr([Xq5(1:mid-1,:);  Xq5(mid+1:end,:)]);
Yq5 =  fliplr([Yq5(1:mid-1,:);  Yq5(mid+1:end,:)]);

cut11 = cut1;
while abs(x(cut11)-Xq4(end,1)) > 1e-9
  cut11++;
endwhile

cut12 = cut11+1;
while abs(x(cut12)-Xq4(1,1)) > 1e-9
  cut12++;
endwhile
clear ans mid


### First Zone of Triangles ###

% reduces the number of cells at the LE down to required at maximum thickness
Xt0 = Xq1(:,end);
Yt0 = Yq1(:,end);

cut3 = cut1;
cut4 = cut2;
N = rows(Xq1);
while N > nLE+1
  cut3++;
  cut4--;
  N = max(N-2,nLE+1);
  Xt0 = [Xt0; interp1([1 N], [x(cut4) x(cut3)], 1:N, "linear")'];
  Yt0 = [Yt0; interp1([1 N], [y(cut4) y(cut3)], 1:N, "linear")'];
endwhile
clear ans N

tri0 = Triangulate(Xt0,Yt0);

### Third Zone of Triangles ###

% reduces the number of cells at the TE to avoid stretched quads
Xt2 = Xq4(:,1);
Yt2 = Yq4(:,1);

cut9  = cut11;
cut10 = cut12;
N = rows(Xq4);
while N > nTAR+1
  cut9--;
  cut10++;
  N = max(N-2,nTAR+1);
  Xt2 = [interp1([1 N], [x(cut10) x(cut9)], 1:N, "linear")'; Xt2];
  Yt2 = [interp1([1 N], [y(cut10) y(cut9)], 1:N, "linear")'; Yt2];
endwhile
clear ans N

tri2 = Triangulate(Xt2,Yt2);

### Third Zone of Quads and Second of Triangles ###

% find the location where the cell size gets to the cell size
% at the maximum thickness points
thresh = data.t*0.01*data.chord/nLE;

cut7 = cut9;  cut8 = cut10;
while true
  cut7--;  cut8++;
  if norm([x(cut7)-x(cut8)  y(cut7)-y(cut8)])/nTAR > thresh
    cut7++;  cut8--;
    break;
  endif
endwhile

clear thresh ans

Xt1 = []; Yt1 = [];

cut5 = cut7+1;  cut6 = cut8-1;
for i=0:nLE-nTAR
  cut5--; cut6++;
  Xt1 = [interp1([1 nTAR+1+i], [x(cut6) x(cut5)], 1:(nTAR+1+i),"linear")'; Xt1];
  Yt1 = [interp1([1 nTAR+1+i], [y(cut6) y(cut5)], 1:(nTAR+1+i),"linear")'; Yt1];
endfor
clear i ans

tri1 = Triangulate(Xt1,Yt1);

% now we can connect the last 2 triangle zones with quads
xl = fliplr(x(cut10:cut8));
yl = fliplr(y(cut10:cut8));

xu = x(cut7:cut9);
yu = y(cut7:cut9);

w = (0:nTAR)'/nTAR;
Xq3 = w*xu + (1-w)*xl;
Yq3 = w*yu + (1-w)*yl;


### Middle (main quad block) ###

xl = fliplr(x(cut6:cut4));
yl = fliplr(y(cut6:cut4));

xu = x(cut3:cut5);
yu = y(cut3:cut5);

w = (0:nLE)'/nLE;
Xq2 = w*xu + (1-w)*xl;
Yq2 = w*yu + (1-w)*yl;

clear xl yl xu yu w


### Create the tip cap ###

% the part of the cap inside the boundary layer is created by morphing the 2D
% cap generated so far, Z is constant for each section. for the second part
% we use interpolation.

% for morphing we use RBF, each outer point is an input point, we do not need
% connectivity information for RBF so for now we do not need to care about "BCs"
% we can either recompute the interpolation matrix based on the previous layer
% or keep it constant, an influence radius of one thickness should be enough

% interpolation kernel
radius = data.t*0.01*data.chord;
xcp = x(1:end-1);
ycp = y(1:end-1);
dist = sqrt((xcp-xcp').^2+(ycp-ycp').^2)/radius;
M = max(1-dist,0).^4.*(1+4*dist);

%distance matrices
xi = reshape(Xq0,numel(Xq0),1);
yi = reshape(Yq0,numel(Yq0),1);
dist = sqrt((xcp-xi).^2+(ycp-yi).^2)/radius;
Aq0 = max(1-dist,0).^4.*(1+4*dist);

xi = reshape(Xq1,numel(Xq1),1);
yi = reshape(Yq1,numel(Yq1),1);
dist = sqrt((xcp-xi).^2+(ycp-yi).^2)/radius;
Aq1 = max(1-dist,0).^4.*(1+4*dist);

xi = reshape(Xq2,numel(Xq2),1);
yi = reshape(Yq2,numel(Yq2),1);
dist = sqrt((xcp-xi).^2+(ycp-yi).^2)/radius;
Aq2 = max(1-dist,0).^4.*(1+4*dist);

xi = reshape(Xq3,numel(Xq3),1);
yi = reshape(Yq3,numel(Yq3),1);
dist = sqrt((xcp-xi).^2+(ycp-yi).^2)/radius;
Aq3 = max(1-dist,0).^4.*(1+4*dist);

xi = reshape(Xq4,numel(Xq4),1);
yi = reshape(Yq4,numel(Yq4),1);
dist = sqrt((xcp-xi).^2+(ycp-yi).^2)/radius;
Aq4 = max(1-dist,0).^4.*(1+4*dist);

xi = reshape(Xq5,numel(Xq5),1);
yi = reshape(Yq5,numel(Yq5),1);
dist = sqrt((xcp-xi).^2+(ycp-yi).^2)/radius;
Aq5 = max(1-dist,0).^4.*(1+4*dist);

dist = sqrt((xcp-Xt0).^2+(ycp-Yt0).^2)/radius;
At0 = max(1-dist,0).^4.*(1+4*dist);

dist = sqrt((xcp-Xt1).^2+(ycp-Yt1).^2)/radius;
At1 = max(1-dist,0).^4.*(1+4*dist);

dist = sqrt((xcp-Xt2).^2+(ycp-Yt2).^2)/radius;
At2 = max(1-dist,0).^4.*(1+4*dist);

clear xi yi dist radius

Zq0 = ones([size(Xq0) BL])*Lz;
Zq1 = ones([size(Xq1) BL])*Lz;
Zq2 = ones([size(Xq2) BL])*Lz;
Zq3 = ones([size(Xq3) BL])*Lz;
Zq4 = ones([size(Xq4) BL])*Lz;
Zq5 = ones([size(Xq5) BL])*Lz;
Zt1 = ones(numel(Xt1),BL)*Lz;
Zt0 = ones(numel(Xt0),BL)*Lz;
Zt2 = ones(numel(Xt2),BL)*Lz;

for i=2:BL
  % control point displacement
  ucp = Xcyl(i,2:end-1,end)' - xcp';
  vcp = Ycyl(i,2:end-1,end)' - ycp';
  wcp = Zcyl(i,2,end) - Lz;
  
  cu = M\ucp;
  cv = M\vcp;
  
  % zone by zone interpolate
  [m, n, ~] = size(Xq0);
  Xq0(:,:,i) = Xq0(:,:,1) + reshape(Aq0*cu, m,n);
  Yq0(:,:,i) = Yq0(:,:,1) + reshape(Aq0*cv, m,n);
  Zq0(:,:,i) = Zq0(:,:,1) + wcp;
  
  [m, n, ~] = size(Xq1);
  Xq1(:,:,i) = Xq1(:,:,1) + reshape(Aq1*cu, m,n);
  Yq1(:,:,i) = Yq1(:,:,1) + reshape(Aq1*cv, m,n);
  Zq1(:,:,i) = Zq1(:,:,1) + wcp;
  
  [m, n, ~] = size(Xq2);
  Xq2(:,:,i) = Xq2(:,:,1) + reshape(Aq2*cu, m,n);
  Yq2(:,:,i) = Yq2(:,:,1) + reshape(Aq2*cv, m,n);
  Zq2(:,:,i) = Zq2(:,:,1) + wcp;
  
  [m, n, ~] = size(Xq3);
  Xq3(:,:,i) = Xq3(:,:,1) + reshape(Aq3*cu, m,n);
  Yq3(:,:,i) = Yq3(:,:,1) + reshape(Aq3*cv, m,n);
  Zq3(:,:,i) = Zq3(:,:,1) + wcp;
  
  [m, n, ~] = size(Xq4);
  Xq4(:,:,i) = Xq4(:,:,1) + reshape(Aq4*cu, m,n);
  Yq4(:,:,i) = Yq4(:,:,1) + reshape(Aq4*cv, m,n);
  Zq4(:,:,i) = Zq4(:,:,1) + wcp;
  
  [m, n, ~] = size(Xq5);
  Xq5(:,:,i) = Xq5(:,:,1) + reshape(Aq5*cu, m,n);
  Yq5(:,:,i) = Yq5(:,:,1) + reshape(Aq5*cv, m,n);
  Zq5(:,:,i) = Zq5(:,:,1) + wcp;
  
  Xt0(:,i) = Xt0(:,1) + At0*cu;
  Yt0(:,i) = Yt0(:,1) + At0*cv;
  Zt0(:,i) = Zt0(:,1) + wcp;
  
  Xt1(:,i) = Xt1(:,1) + At1*cu;
  Yt1(:,i) = Yt1(:,1) + At1*cv;
  Zt1(:,i) = Zt1(:,1) + wcp;
  
  Xt2(:,i) = Xt2(:,1) + At2*cu;
  Yt2(:,i) = Yt2(:,1) + At2*cv;
  Zt2(:,i) = Zt2(:,1) + wcp;
endfor

clear i m n Aq0 Aq1 Aq2 Aq3 Aq4 Aq5 At0 At1 At2 M cu cv ucp vcp wcp xcp ycp

% map the last section of the cap to the surface of the sphere and continue the
% extrusion using the spacing determined for the last section of the annulus.

for i=2:numel(R3)
  Xq0(:,:,end+1) = 0;  Yq0(:,:,end+1) = 0;  Zq0(:,:,end+1) = 0;
  Xq1(:,:,end+1) = 0;  Yq1(:,:,end+1) = 0;  Zq1(:,:,end+1) = 0;
  Xq2(:,:,end+1) = 0;  Yq2(:,:,end+1) = 0;  Zq2(:,:,end+1) = 0;
  Xq3(:,:,end+1) = 0;  Yq3(:,:,end+1) = 0;  Zq3(:,:,end+1) = 0;
  Xq4(:,:,end+1) = 0;  Yq4(:,:,end+1) = 0;  Zq4(:,:,end+1) = 0;
  Xq5(:,:,end+1) = 0;  Yq5(:,:,end+1) = 0;  Zq5(:,:,end+1) = 0;
  Xt0(: , end+1) = 0;  Yt0(: , end+1) = 0;  Zt0(: , end+1) = 0;
  Xt1(: , end+1) = 0;  Yt1(: , end+1) = 0;  Zt1(: , end+1) = 0;
  Xt2(: , end+1) = 0;  Yt2(: , end+1) = 0;  Zt2(: , end+1) = 0;
endfor

% use the same scaled projection used for the annulus
phi = atan2(Yq0(:,:,BL),Xq0(:,:,BL));
theta = atan(sqrt(Yq0(:,:,BL).^2+Xq0(:,:,BL).^2)/focus)*data.projFactor;
Xq0(:,:,end) = Rsph * sin(theta).* cos(phi);
Yq0(:,:,end) = Rsph * sin(theta).* sin(phi);
Zq0(:,:,end) = Rsph * cos(theta);

phi = atan2(Yq1(:,:,BL),Xq1(:,:,BL));
theta = atan(sqrt(Yq1(:,:,BL).^2+Xq1(:,:,BL).^2)/focus)*data.projFactor;
Xq1(:,:,end) = Rsph * sin(theta).* cos(phi);
Yq1(:,:,end) = Rsph * sin(theta).* sin(phi);
Zq1(:,:,end) = Rsph * cos(theta);

phi = atan2(Yq2(:,:,BL),Xq2(:,:,BL));
theta = atan(sqrt(Yq2(:,:,BL).^2+Xq2(:,:,BL).^2)/focus)*data.projFactor;
Xq2(:,:,end) = Rsph * sin(theta).* cos(phi);
Yq2(:,:,end) = Rsph * sin(theta).* sin(phi);
Zq2(:,:,end) = Rsph * cos(theta);

phi = atan2(Yq3(:,:,BL),Xq3(:,:,BL));
theta = atan(sqrt(Yq3(:,:,BL).^2+Xq3(:,:,BL).^2)/focus)*data.projFactor;
Xq3(:,:,end) = Rsph * sin(theta).* cos(phi);
Yq3(:,:,end) = Rsph * sin(theta).* sin(phi);
Zq3(:,:,end) = Rsph * cos(theta);

phi = atan2(Yq4(:,:,BL),Xq4(:,:,BL));
theta = atan(sqrt(Yq4(:,:,BL).^2+Xq4(:,:,BL).^2)/focus)*data.projFactor;
Xq4(:,:,end) = Rsph * sin(theta).* cos(phi);
Yq4(:,:,end) = Rsph * sin(theta).* sin(phi);
Zq4(:,:,end) = Rsph * cos(theta);

phi = atan2(Yq5(:,:,BL),Xq5(:,:,BL));
theta = atan(sqrt(Yq5(:,:,BL).^2+Xq5(:,:,BL).^2)/focus)*data.projFactor;
Xq5(:,:,end) = Rsph * sin(theta).* cos(phi);
Yq5(:,:,end) = Rsph * sin(theta).* sin(phi);
Zq5(:,:,end) = Rsph * cos(theta);

phi = atan2(Yt0(:,BL),Xt0(:,BL));
theta = atan(sqrt(Yt0(:,BL).^2+Xt0(:,BL).^2)/focus)*data.projFactor;
Xt0(:,end) = Rsph * sin(theta).* cos(phi);
Yt0(:,end) = Rsph * sin(theta).* sin(phi);
Zt0(:,end) = Rsph * cos(theta);

phi = atan2(Yt1(:,BL),Xt1(:,BL));
theta = atan(sqrt(Yt1(:,BL).^2+Xt1(:,BL).^2)/focus)*data.projFactor;
Xt1(:,end) = Rsph * sin(theta).* cos(phi);
Yt1(:,end) = Rsph * sin(theta).* sin(phi);
Zt1(:,end) = Rsph * cos(theta);

phi = atan2(Yt2(:,BL),Xt2(:,BL));
theta = atan(sqrt(Yt2(:,BL).^2+Xt2(:,BL).^2)/focus)*data.projFactor;
Xt2(:,end) = Rsph * sin(theta).* cos(phi);
Yt2(:,end) = Rsph * sin(theta).* sin(phi);
Zt2(:,end) = Rsph * cos(theta);

for i=1:numel(R3)-2
  w = R3(i+1);
  
  Xq0(:,:,BL+i) = (1-w)*Xq0(:,:,BL) + w*Xq0(:,:,end);
  Yq0(:,:,BL+i) = (1-w)*Yq0(:,:,BL) + w*Yq0(:,:,end);
  Zq0(:,:,BL+i) = (1-w)*Zq0(:,:,BL) + w*Zq0(:,:,end);
  
  Xq1(:,:,BL+i) = (1-w)*Xq1(:,:,BL) + w*Xq1(:,:,end);
  Yq1(:,:,BL+i) = (1-w)*Yq1(:,:,BL) + w*Yq1(:,:,end);
  Zq1(:,:,BL+i) = (1-w)*Zq1(:,:,BL) + w*Zq1(:,:,end);
  
  Xq2(:,:,BL+i) = (1-w)*Xq2(:,:,BL) + w*Xq2(:,:,end);
  Yq2(:,:,BL+i) = (1-w)*Yq2(:,:,BL) + w*Yq2(:,:,end);
  Zq2(:,:,BL+i) = (1-w)*Zq2(:,:,BL) + w*Zq2(:,:,end);
  
  Xq3(:,:,BL+i) = (1-w)*Xq3(:,:,BL) + w*Xq3(:,:,end);
  Yq3(:,:,BL+i) = (1-w)*Yq3(:,:,BL) + w*Yq3(:,:,end);
  Zq3(:,:,BL+i) = (1-w)*Zq3(:,:,BL) + w*Zq3(:,:,end);
  
  Xq4(:,:,BL+i) = (1-w)*Xq4(:,:,BL) + w*Xq4(:,:,end);
  Yq4(:,:,BL+i) = (1-w)*Yq4(:,:,BL) + w*Yq4(:,:,end);
  Zq4(:,:,BL+i) = (1-w)*Zq4(:,:,BL) + w*Zq4(:,:,end);
  
  Xq5(:,:,BL+i) = (1-w)*Xq5(:,:,BL) + w*Xq5(:,:,end);
  Yq5(:,:,BL+i) = (1-w)*Yq5(:,:,BL) + w*Yq5(:,:,end);
  Zq5(:,:,BL+i) = (1-w)*Zq5(:,:,BL) + w*Zq5(:,:,end);
  
  Xt0(: , BL+i) = (1-w)*Xt0(: , BL) + w*Xt0(: , end);
  Yt0(: , BL+i) = (1-w)*Yt0(: , BL) + w*Yt0(: , end);
  Zt0(: , BL+i) = (1-w)*Zt0(: , BL) + w*Zt0(: , end);
  
  Xt1(: , BL+i) = (1-w)*Xt1(: , BL) + w*Xt1(: , end);
  Yt1(: , BL+i) = (1-w)*Yt1(: , BL) + w*Yt1(: , end);
  Zt1(: , BL+i) = (1-w)*Zt1(: , BL) + w*Zt1(: , end);
  
  Xt2(: , BL+i) = (1-w)*Xt2(: , BL) + w*Xt2(: , end);
  Yt2(: , BL+i) = (1-w)*Yt2(: , BL) + w*Yt2(: , end);
  Zt2(: , BL+i) = (1-w)*Zt2(: , BL) + w*Zt2(: , end);
endfor

clear i w phi theta R3 focus

% smooth corners of the leading and trailing edge blocks
w = 0.7;

U = (Xq0(end-1,end,:) + Xq0(end,end-1,:) + ...
     Xq1(end+1-columns(Xq0),2,:))/3 - Xq0(end,end,:);
V = (Yq0(end-1,end,:) + Yq0(end,end-1,:) + ...
     Yq1(end+1-columns(Xq0),2,:))/3 - Yq0(end,end,:);
W = (Zq0(end-1,end,:) + Zq0(end,end-1,:) + ...
     Zq1(end+1-columns(Xq0),2,:))/3 - Zq0(end,end,:);

Xq0(end,end,:) += w*U;  Xq1(end+1-columns(Xq0),1,:) += w*U;
Yq0(end,end,:) += w*V;  Yq1(end+1-columns(Xq0),1,:) += w*V;
Zq0(end,end,:) += w*W;  Zq1(end+1-columns(Xq0),1,:) += w*W;


U = (Xq0(2,end,:) + Xq0(1,end-1,:) + Xq1(columns(Xq0),2,:))/3 - Xq0(1,end,:);
V = (Yq0(2,end,:) + Yq0(1,end-1,:) + Yq1(columns(Xq0),2,:))/3 - Yq0(1,end,:);
W = (Zq0(2,end,:) + Zq0(1,end-1,:) + Zq1(columns(Xq0),2,:))/3 - Zq0(1,end,:);

Xq0(1,end,:) += w*U;  Xq1(columns(Xq0),1,:) += w*U;
Yq0(1,end,:) += w*V;  Yq1(columns(Xq0),1,:) += w*V;
Zq0(1,end,:) += w*W;  Zq1(columns(Xq0),1,:) += w*W;


U = (Xq5(end-1,1,:) + Xq5(end,2,:) + ...
     Xq4(end+1-columns(Xq5),end-1,:))/3 - Xq5(end,1,:);
V = (Yq5(end-1,1,:) + Yq5(end,2,:) + ...
     Yq4(end+1-columns(Xq5),end-1,:))/3 - Yq5(end,1,:);
W = (Zq5(end-1,1,:) + Zq5(end,2,:) + ...
     Zq4(end+1-columns(Xq5),end-1,:))/3 - Zq5(end,1,:);

Xq5(end,1,:) += w*U;  Xq4(end+1-columns(Xq5),end,:) += w*U;
Yq5(end,1,:) += w*V;  Yq4(end+1-columns(Xq5),end,:) += w*V;
Zq5(end,1,:) += w*W;  Zq4(end+1-columns(Xq5),end,:) += w*W;


U = (Xq5(2,1,:) + Xq5(1,2,:) + Xq4(columns(Xq5),end-1,:))/3 - Xq5(1,1,:);
V = (Yq5(2,1,:) + Yq5(1,2,:) + Yq4(columns(Xq5),end-1,:))/3 - Yq5(1,1,:);
W = (Zq5(2,1,:) + Zq5(1,2,:) + Zq4(columns(Xq5),end-1,:))/3 - Zq5(1,1,:);

Xq5(1,1,:) += w*U;  Xq4(columns(Xq5),end,:) += w*U;
Yq5(1,1,:) += w*V;  Yq4(columns(Xq5),end,:) += w*V;
Zq5(1,1,:) += w*W;  Zq4(columns(Xq5),end,:) += w*W;

clear U V W Rsph w x y

% save block sections to create solid domain
% before any deformation is applied to the fluid domain
if data.fsi
  section.xq0 = Xq0(:,:,1);  section.yq0 = Yq0(:,:,1);
  section.xq1 = Xq1(:,:,1);  section.yq1 = Yq1(:,:,1);
  section.xq2 = Xq2(:,:,1);  section.yq2 = Yq2(:,:,1);
  section.xq3 = Xq3(:,:,1);  section.yq3 = Yq3(:,:,1);
  section.xq4 = Xq4(:,:,1);  section.yq4 = Yq4(:,:,1);
  section.xq5 = Xq5(:,:,1);  section.yq5 = Yq5(:,:,1);

  section.xt0 = Xt0(:,1);  section.yt0 = Yt0(:,1);
  section.xt1 = Xt1(:,1);  section.yt1 = Yt1(:,1);
  section.xt2 = Xt2(:,1);  section.yt2 = Yt2(:,1);
endif
