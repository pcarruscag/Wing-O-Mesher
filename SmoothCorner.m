%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
% Smooth the wing tip corner
if(!data.smoothTip) return; endif

% move points along the 45 degree line of the tip cap to reduce skewness
function [xnew, ynew, znew] = deform(xold, yold, zold, xref, yref, zref, Xnew)
  % map old to reference points
  [m,n,p] = size(xold);
  s = m*n*p;
  xold = reshape(xold,s,1);
  yold = reshape(yold,s,1);
  zold = reshape(zold,s,1);
  tol = 1e-7;
  map = (abs(xold-xref)<tol & abs(yold-yref)<tol &...
         abs(zold-zref)<tol) * (1:numel(xref))';
  hasMatch = map!=0;
  map = map(hasMatch);
  xold(hasMatch) = Xnew(1,map)';
  yold(hasMatch) = Xnew(2,map)';
  zold(hasMatch) = Xnew(3,map)';
  xnew = reshape(xold,m,n,p);
  ynew = reshape(yold,m,n,p);
  znew = reshape(zold,m,n,p);
endfunction

for i=2:BL
  w = 1-((i-2)/(BL-2))^4;
  % initial distance
  di = norm([X(i,:)-X(1,:); Y(i,:)-Y(1,:)],"cols");
  % current distance
  dir = [Xcyl(i,:,end)-Xcyl(1,:,end); Ycyl(i,:,end)-Ycyl(1,:,end);...
         Zcyl(i,:,end)-Zcyl(1,:,end)];
  dc = norm(dir,"cols");
  % normalized direction
  dir = dir ./ dc;
  % target distance
  dt = w*di + (1-w)*dc;
  % new position
  pos = dir .* dt + [Xcyl(1,:,end); Ycyl(1,:,end); Zcyl(1,:,end)];
  
  % apply new position
  xref = Xcyl(i,2:end-1,end);
  yref = Ycyl(i,2:end-1,end);
  zref = Zcyl(i,2:end-1,end);
  Xcyl(i,:,end) = pos(1,:);
  Ycyl(i,:,end) = pos(2,:);
  Zcyl(i,:,end) = pos(3,:);
  pos = pos(:,2:end-1);
  [Xq0(:,:,i),Yq0(:,:,i),Zq0(:,:,i)]=deform(Xq0(:,:,i),Yq0(:,:,i),Zq0(:,:,i),...
                                            xref, yref, zref ,pos);
  [Xq1(:,:,i),Yq1(:,:,i),Zq1(:,:,i)]=deform(Xq1(:,:,i),Yq1(:,:,i),Zq1(:,:,i),...
                                            xref, yref, zref ,pos);
  [Xq2(:,:,i),Yq2(:,:,i),Zq2(:,:,i)]=deform(Xq2(:,:,i),Yq2(:,:,i),Zq2(:,:,i),...
                                            xref, yref, zref ,pos);
  [Xq3(:,:,i),Yq3(:,:,i),Zq3(:,:,i)]=deform(Xq3(:,:,i),Yq3(:,:,i),Zq3(:,:,i),...
                                            xref, yref, zref ,pos);
  [Xq4(:,:,i),Yq4(:,:,i),Zq4(:,:,i)]=deform(Xq4(:,:,i),Yq4(:,:,i),Zq4(:,:,i),...
                                            xref, yref, zref ,pos);
  [Xq5(:,:,i),Yq5(:,:,i),Zq5(:,:,i)]=deform(Xq5(:,:,i),Yq5(:,:,i),Zq5(:,:,i),...
                                            xref, yref, zref ,pos);
  [Xt0(:,i) , Yt0(:,i) , Zt0(:,i)] = deform(Xt0(:,i), Yt0(:,i), Zt0(:,i),...
                                            xref, yref, zref ,pos);
  [Xt1(:,i) , Yt1(:,i) , Zt1(:,i)] = deform(Xt1(:,i), Yt1(:,i), Zt1(:,i),...
                                            xref, yref, zref ,pos);
  [Xt2(:,i) , Yt2(:,i) , Zt2(:,i)] = deform(Xt2(:,i), Yt2(:,i), Zt2(:,i),...
                                            xref, yref, zref ,pos);
endfor

Xan(1,:,end) = Xcyl(BL,:,end);
Yan(1,:,end) = Ycyl(BL,:,end);
Zan(1,:,end) = Zcyl(BL,:,end);

clear i w di dir dc dt pos xref yref zref

% consider this many sections close to the tip when building an RBF
endOffset = 5;
moveFactor = [0.35 1.25 0.15 0.3];

[Xwing,Ywing,Zwing,Nwing] = GetWingSurfPoints(endOffset, Xcyl,Ycyl,Zcyl,...
        Xq0,Yq0,Zq0, Xq1,Yq1,Zq1, Xt0,Yt0,Zt0, Xq2,Yq2,Zq2, Xt1,Yt1,Zt1,...
        Xq3,Yq3,Zq3, Xt2,Yt2,Zt2, Xq4,Yq4,Zq4, Xq5,Yq5,Zq5);

%build list of inner and outer points for the last section
% coordinates of surface points, these were not included in this order in "wing"
Xsurf = squeeze(Xcyl(1,2:end-1,end));
Ysurf = squeeze(Ycyl(1,2:end-1,end));

% coordinates of tip section, including "surf" in unknown order
Xinner = Xwing(Nwing+1:end);
Yinner = Ywing(Nwing+1:end);

% match "surf" and "tip"
tol = 1e-7;
tmp = abs(Xinner-Xsurf)<tol & abs(Yinner-Ysurf)<tol;
assert(sum(sum(tmp))==numel(Xsurf),"could not find outer layer");
outerMask = any(tmp,2);

% split "inner/outer" in the order they appear in "wing", map "outer" to "wing"
mapOut = (1:numel(Xinner))(outerMask)+Nwing;
mapIn  = (1:numel(Xinner))(!outerMask)+Nwing;
Xouter = Xinner(outerMask);  Xinner = Xinner(!outerMask);
Youter = Yinner(outerMask);  Yinner = Yinner(!outerMask);

% for each "surf" point find the closest "inner" point
tmp = (Xinner-Xsurf).^2 + (Yinner-Ysurf).^2;
[tmp,idx] = min(tmp);
d = sqrt(tmp);
% smooth distance for z displacement (would not work with "outer")
tmp = [d d d];
for i=1:20
  tmp = [0 0 (tmp(1:end-4)+tmp(2:end-3)+tmp(3:end-2)+...
         tmp(4:end-1)+tmp(5:end))*0.2 0 0];
endfor
d = tmp(numel(d)+1:2*numel(d));

% permute results as map is for "outer" not "surf"
tmp = abs(Xouter-Xsurf)<tol & abs(Youter-Ysurf)<tol;
idx *= tmp';
d *= tmp';

clear tol i tmp Xsurf Ysurf outerMask

% displacements to apply to outer points of last section to smooth corner
move = moveFactor(1)*data.taper;
usmooth =  move*(Xinner(idx)-Xouter);
vsmooth =  move*(Yinner(idx)-Youter);
wsmooth = -min(move*d'*moveFactor(2), 0.5*(Zwing(Nwing+1)-Zwing(Nwing)));

u = zeros(size(Xwing)); v = u; w = u;
u(mapOut) += usmooth;  u(mapIn(idx)) += moveFactor(3)*usmooth;
v(mapOut) += vsmooth;  v(mapIn(idx)) += moveFactor(3)*vsmooth;
w(mapOut) += wsmooth;  w(mapIn(idx)) += moveFactor(4)*wsmooth;

clear Xinner Yinner Xouter Youter idx d mapOut
clear usmooth vsmooth wsmooth move Nwing mapIn

% Create RBF kernel to apply displacements to the inner mesh points.
radius = 3*max(abs(u));
tmp = SolveRBF(Xwing/radius, Ywing/radius, Zwing/radius, [u,v,w]);

% store state needed to apply deformation
smoother.radius = radius;
smoother.xcp = Xwing;  smoother.ucp = tmp(:,1);
smoother.ycp = Ywing;  smoother.vcp = tmp(:,2);
smoother.zcp = Zwing;  smoother.wcp = tmp(:,3);

clear u v w dist M tmp Xwing Ywing Zwing

% Apply displacements
fun = defXYZ; % rounding the corner affects all coordinates

% use structure of the grid to avoid computing distances to all points
R = sqrt((X-X(1,:)).^2+(Y-Y(1,:)).^2);
idx_cyl = sum(min(R,[],2) < radius)+1;

mid = ceil(rows(Xq0)/2);
R = sqrt((Xq0(mid,1,:)-Xq0(mid,1,1)).^2+...
         (Yq0(mid,1,:)-Yq0(mid,1,1)).^2+...
         (Zq0(mid,1,:)-Zq0(mid,1,1)).^2);
idx_tip = sum(squeeze(R) < radius)+1;

% this deformation is very local, annulus should not be affected
assert(idx_tip-data.BL < 0, "deformation is too large");

clear R mid

% cylinder
[m,n,p] = size(Xcyl(1:idx_cyl,:,end-endOffset:end));
x = reshape(Xcyl(1:idx_cyl,:,end-endOffset:end),m*n*p,1);
y = reshape(Ycyl(1:idx_cyl,:,end-endOffset:end),m*n*p,1);
z = reshape(Zcyl(1:idx_cyl,:,end-endOffset:end),m*n*p,1);
[x,y,z] = fun(x,y,z,smoother);
Xcyl(1:idx_cyl,:,end-endOffset:end) = reshape(x,m,n,p)*smoother.radius;
Ycyl(1:idx_cyl,:,end-endOffset:end) = reshape(y,m,n,p)*smoother.radius;
Zcyl(1:idx_cyl,:,end-endOffset:end) = reshape(z,m,n,p)*smoother.radius;

% tip extrusion
% hex blocks
[Xq0,Yq0,Zq0] = defHexBlock(Xq0,Yq0,Zq0,fun,smoother,idx_tip);
[Xq1,Yq1,Zq1] = defHexBlock(Xq1,Yq1,Zq1,fun,smoother,idx_tip);
[Xq2,Yq2,Zq2] = defHexBlock(Xq2,Yq2,Zq2,fun,smoother,idx_tip);
[Xq3,Yq3,Zq3] = defHexBlock(Xq3,Yq3,Zq3,fun,smoother,idx_tip);
[Xq4,Yq4,Zq4] = defHexBlock(Xq4,Yq4,Zq4,fun,smoother,idx_tip);
[Xq5,Yq5,Zq5] = defHexBlock(Xq5,Yq5,Zq5,fun,smoother,idx_tip);

%prism blocks
[Xt0,Yt0,Zt0] = defPrismBlock(Xt0,Yt0,Zt0,fun,smoother,idx_tip);
[Xt1,Yt1,Zt1] = defPrismBlock(Xt1,Yt1,Zt1,fun,smoother,idx_tip);
[Xt2,Yt2,Zt2] = defPrismBlock(Xt2,Yt2,Zt2,fun,smoother,idx_tip);

clear m n p x y z ans idx_cyl idx_tip endOffset moveFactor radius chunks fun
if(!data.fsi) clear("smoother","defXYZ"); endif
