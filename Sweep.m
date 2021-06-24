%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
% Apply sweep angle (of the 0.25c line) and taper to the wing via RBF deform.
if(!data.sweep && data.taper==1 && !data.twist && !data.extDef) return; endif

% control points before motion
c = data.chord;
t = data.t*0.01*c;
CP0 = 0.5*[-c  -t;...
            c  -t;...
           -c   t;...
            c   t];
dc = 0.5*(1-data.taper)*c;
dt = 0.5*(1-data.taper*data.thinning)*t;
dCP = [ dc   dt;...
       -dc   dt;...
        dc  -dt;...
       -dc  -dt];

% control points after taper
CP1 = CP0 + dCP;

% apply twist
rotMat = [cos(data.twist) -sin(data.twist);...
          sin(data.twist)  cos(data.twist)];
CP1 = CP1*rotMat;

% translation of 0.25c point for required sweep
dx = tan(data.sweep*pi/180)*Lz;
% correct for what is already induced by taper
dx -= 0.25*(1-data.taper)*data.chord;

CP1(:,1) += dx;
dCP = CP1 - CP0;

clear dc dt dx rotMat

% dCP are x and y displacements of the control points (there is no z motion).
% We use the basis ("shape") functions of our trivial box to map them to the
% surface of the wing.

[Xwing,Ywing,Zwing] = GetWingSurfPoints(1000, Xcyl,Ycyl,Zcyl, Xq0,Yq0,Zq0,...
                       Xq1,Yq1,Zq1, Xt0,Yt0,Zt0, Xq2,Yq2,Zq2, Xt1,Yt1,Zt1,...
                       Xq3,Yq3,Zq3, Xt2,Yt2,Zt2, Xq4,Yq4,Zq4, Xq5,Yq5,Zq5);
assert(numel(Xwing)**2/2^28 < RAM_LIMIT,'Stopping to not go above RAM limit.');
Xi = 2*Xwing/c;
Eta = 2*Ywing/t;

P1 = dCP(1,1);  P2 = dCP(2,1);  P3 = dCP(3,1);  P4 = dCP(4,1);

u = P1*(1-Xi).*(1-Eta) + P2*(1+Xi).*(1-Eta)+...
    P3*(1-Xi).*(1+Eta) + P4*(1+Xi).*(1+Eta);
u = u.*Zwing*0.25/Lz;

P1 = dCP(1,2);  P2 = dCP(2,2);  P3 = dCP(3,2);  P4 = dCP(4,2);

v = P1*(1-Xi).*(1-Eta) + P2*(1+Xi).*(1-Eta)+...
    P3*(1-Xi).*(1+Eta) + P4*(1+Xi).*(1+Eta);
v = v.*Zwing*0.25/Lz;

clear c t Xi Eta P1 P2 P3 P4

% share u deformation between root and tip to avoid excessive skew
ave_du = 0.5*(min(min(dCP(:,1)),0)+max(max(dCP(:,1)),0));
u -= ave_du;
CP0(:,1) -= ave_du;
CP1(:,1) -= ave_du;

% apply external deformation on top of the simple linear taper and twist
if data.extDef
  % Compute deformation based on reference and deformed points,
  % the reference should match the "current" shape of the wing.

  ref = csvread("external_deformation/reference.csv")(:,[2:end]);
  assert(numel(Xwing) == rows(ref),'Sizes do not match.');
  def = csvread("external_deformation/deformed.csv")(:,[2,4]) - ref(:,[1,3]);
  % swap y/z
  tmp = ref(:,2);
  ref(:,2) = ref(:,3);
  ref(:,3) = tmp;
  clear tmp;

  % compute "current" shape and map it to the external reference
  curr = [Xwing+u, Ywing+v, Zwing];
  map = zeros(rows(ref),1);
  chunk = 1024;
  cursor = 1;
  while cursor <= rows(ref)
    idx = cursor : min(cursor+chunk, rows(ref));
    [~,map(idx)] = min(abs(curr(idx,1)-ref(:,1)')+...
                       abs(curr(idx,2)-ref(:,2)')+...
                       abs(curr(idx,3)-ref(:,3)'),[],2);
    cursor += chunk;
  endwhile
  assert(numel(map) == numel(unique(map)),'Geometries do not match.');

  % update displacements
  u += def(map,1);
  v += def(map,2);
  clear ref def curr map chunk cursor idx

  % re-"center" deformation to reduce skew (give more importance to tip)
  ave_du = 0.5*(max(u) + min(u));
  u -= ave_du;
  CP0(:,1) -= ave_du;
  CP1(:,1) -= ave_du;
endif

% Create RBF kernel to apply displacements to the inner mesh points.
radius = 4*max(max(abs(u)), max(abs(v)));
tmp = SolveRBF(Xwing/radius, Ywing/radius, Zwing/radius, [u,v]);

% store state needed to apply deformation
sweeper.radius = radius;
sweeper.xcp = Xwing;  sweeper.ucp = tmp(:,1);
sweeper.ycp = Ywing;  sweeper.vcp = tmp(:,2);
sweeper.zcp = Zwing;

clear u v dCP dist M tmp Xwing Ywing Zwing fid ave_du Lz

% Apply displacements
fun = defXY; % sweeping does not modify z

% use structure of the grid to avoid computing distances to all points
R = sqrt((X-X(1,:)).^2+(Y-Y(1,:)).^2);
idx_cyl = sum(min(R,[],2) < radius)+1;

mid = ceil(rows(Xq0)/2);
R = sqrt((Xq0(mid,1,:)-Xq0(mid,1,1)).^2+...
         (Yq0(mid,1,:)-Yq0(mid,1,1)).^2+...
         (Zq0(mid,1,:)-Zq0(mid,1,1)).^2);
idx_tip = sum(squeeze(R) < radius)+1;

idx_an1 = max(1,idx_tip-data.BL+1);
idx_an2 = size(Xan)(3)-idx_cyl+data.BL;

clear R mid

% cylinder
[m,n,p] = size(Xcyl(1:idx_cyl,:,:));
x = reshape(Xcyl(1:idx_cyl,:,:),m*n*p,1);
y = reshape(Ycyl(1:idx_cyl,:,:),m*n*p,1);
z = reshape(Zcyl(1:idx_cyl,:,:),m*n*p,1);
[x,y,z] = fun(x,y,z,sweeper);
Xcyl(1:idx_cyl,:,:) = reshape(x,m,n,p)*sweeper.radius;
Ycyl(1:idx_cyl,:,:) = reshape(y,m,n,p)*sweeper.radius;

% annulus
[m,n,p] = size(Xan(1:idx_an1,:,idx_an2:end));
x = reshape(Xan(1:idx_an1,:,idx_an2:end),m*n*p,1);
y = reshape(Yan(1:idx_an1,:,idx_an2:end),m*n*p,1);
z = reshape(Zan(1:idx_an1,:,idx_an2:end),m*n*p,1);
[x,y,z] = fun(x,y,z,sweeper);
Xan(1:idx_an1,:,idx_an2:end) = reshape(x,m,n,p)*sweeper.radius;
Yan(1:idx_an1,:,idx_an2:end) = reshape(y,m,n,p)*sweeper.radius;

% tip extrusion
% hex blocks
[Xq0,Yq0] = defHexBlock(Xq0,Yq0,Zq0,fun,sweeper,idx_tip);
[Xq1,Yq1] = defHexBlock(Xq1,Yq1,Zq1,fun,sweeper,idx_tip);
[Xq2,Yq2] = defHexBlock(Xq2,Yq2,Zq2,fun,sweeper,idx_tip);
[Xq3,Yq3] = defHexBlock(Xq3,Yq3,Zq3,fun,sweeper,idx_tip);
[Xq4,Yq4] = defHexBlock(Xq4,Yq4,Zq4,fun,sweeper,idx_tip);
[Xq5,Yq5] = defHexBlock(Xq5,Yq5,Zq5,fun,sweeper,idx_tip);

%prism blocks
[Xt0,Yt0] = defPrismBlock(Xt0,Yt0,Zt0,fun,sweeper,idx_tip);
[Xt1,Yt1] = defPrismBlock(Xt1,Yt1,Zt1,fun,sweeper,idx_tip);
[Xt2,Yt2] = defPrismBlock(Xt2,Yt2,Zt2,fun,sweeper,idx_tip);

clear m n p x y z ans idx_cyl idx_tip idx_an1 idx_an2 fun radius chunks
if(!data.fsi) clear("sweeper","defXY"); endif
