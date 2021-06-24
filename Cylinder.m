%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
% Build a cylinder from the O-Grid.
% Last layer has a 45 degree boundary layer zone so that a tip cap can be put
% over it, that will avoid mapping very small cells onto the half sphere.
% Push the corner of the cylinder inwards to trade skewness of the extensions
% for skewness inside the cylinder,

Nz = data.Nz;
BL = data.BL;

Lz = data.Lz*data.chord;
t = (0:Nz-1)/(Nz-1);
if strcmp(data.tipType,'sin')
  k = (data.tipFactor*data.t*data.chord*0.01/data.nTh1/Lz-t(2))/...
      (t(end-1)-sin(0.5*pi*t(end-1)));
  k = min(max(k,0),1);
  z = Lz*((1-k)*t + k*sin(0.5*pi*t));
  t = Lz*(0.5*t + 0.5*sin(0.5*pi*t))*(1-data.dL);
  clear k
elseif strcmp(data.tipType,'tanh')
  z = Lz*tanh(data.tipFactor*t)/tanh(data.tipFactor);
  t = Lz*tanh(0.75*data.tipFactor*t)/tanh(0.75*data.tipFactor)*(1-data.dL);
  % reduce spacing close to the symmetry plane
  if data.rootFactor < 1.0
    dz = z(2);
    z(2) = dz * data.rootFactor;
    for i=3:numel(z)
      zi = z(i-1) + data.rate*(z(i-1)-z(i-2));
      if zi-z(i-1) > dz
        i -= 1;
        break;
      endif
      z(i) = zi;
    endfor
    nz = Nz-i+1;
    nz = (Lz-z(i))*tanh(data.tipFactor*(0:nz-1)/(nz-1))/tanh(data.tipFactor);
    z(i:end) = z(i) + nz;
    clear dx i zi nz
  endif
else
  z = t*Lz;
  t = z*(1-data.dL);
endif

Xcyl = repmat(X,1,1,Nz);
Ycyl = repmat(Y,1,1,Nz);
Zcyl = zeros(data.M,data.N+1,Nz);
for i=1:data.M
  % away from the wall we want linear spacing
  w = (norm([X(i,2)-X(1,2)  Y(i,2)-Y(1,2)])/...
       norm([X(end,2)-X(1,2)  Y(end,2)-Y(1,2)]))^2;
  Zcyl(i,:,:) = repmat((1-w)*z + w*t, data.N+1, 1);
endfor
clear i t w dz

% offset last layer
for i=2:BL
  offset = mean(norm([X(i,:)-X(1,:); Y(i,:)-Y(1,:)],"cols"));
  Zcyl(i,:,end) += offset;
endfor
Zcyl(BL+1:end,:,end) += offset;
clear offset i

% blend the offset
for i=1:3
  w = (z(end-i)-z(end-i-1))/(z(end-i+1)-z(end-i-1));
  Zcyl(:,:,end-i) = w*Zcyl(:,:,end-i+1)+(1-w)*Zcyl(:,:,end-i-1);
endfor
Zcyl(BL,:,end-3:end-1)=0.5*(Zcyl(BL-1,:,end-3:end-1)+Zcyl(BL+1,:,end-3:end-1));
clear i w z

if data.dR != 0.0
  % taper the cylinder
  % interpolation kernel
  radius = 0.75*data.L*data.chord;
  xcp = X(end,2:end-1);
  ycp = Y(end,2:end-1);
  dist = sqrt((xcp-xcp').^2+(ycp-ycp').^2)/radius;
  M = max(1-dist,0).^4.*(1+4*dist);

  %distance matrix
  xi = reshape(X,numel(X),1);
  yi = reshape(Y,numel(Y),1);
  dist = sqrt((xcp-xi).^2+(ycp-yi).^2)/radius;
  A = max(1-dist,0).^4.*(1+4*dist);

  % direction of deformation
  theta = atan2(ycp,xcp);
  ucp = cos(theta);
  vcp = sin(theta);
  cu = M\ucp';
  cv = M\vcp';
  u = reshape(A*cu,data.M,data.N+1);
  v = reshape(A*cv,data.M,data.N+1);

  % deform
  for i=2:Nz
    w = ((i-1)/(Nz-1))^2*data.dR*data.L*data.chord;
    Xcyl(:,:,i) -= w*u;
    Ycyl(:,:,i) -= w*v;
  endfor

  clear radius xcp ycp dist M xi yi A theta ucp vcp cu cv u v i w
endif

% Extend the cylinder by mapping its outer layer to a half sphere
Rsph = data.Rsph*Lz;

% linear distribution of polar angle
theta = 0.5*pi*(1:-0.5/(Nz-1):0.5);

% infer distribution of azimuthal angle
phi = atan2(Y(end,:), X(end,:));

% radial spacing
% locations based on aspect ratio and expansion
R1 = [norm([X(end-1,1) Y(end-1,1)])  data.L*data.chord];
while (Rsph-R1(end))/(R1(end)-R1(end-1)) > 0.5
  i = numel(R1);
  % constant rate
  dR = data.sphRate*(R1(i)-R1(i-1));
  % limit aspect ratio
  maxAR = [data.farAR data.sphAR](1+(R1(end)>data.targetFar*data.chord));
  dT = 2*pi*R1(i)/(data.N-1);
  dR = min(dR,dT*maxAR);
  % append to R
  R1 = [R1 R1(end)+dR];
endwhile
R1(end) = Rsph;
R1(end-1) = 0.5*(Rsph+R1(end-2));
clear i dR dT

% convert to non-dimensional spacing
R1 = (R1(2:end)-R1(2))/(Rsph-R1(2));

% the spacing rule needs to vary with theta since on the last section we are
% extruding the boundary layer. on the last "z" slice we use a power law.

% initial cell size
ds0 = data.rate^2*mean(norm([X(BL+1,:)-X(BL,:); Y(BL+1,:)-Y(BL,:)],"cols"));

% representative total length we need to cover
L = Rsph-Lz;
r0 = data.sphRate;
for i=1:100
  r = (1-L*(1-r0)/ds0)^(1/(numel(R1)-1));
  if abs(r-r0) < 1e-8
    break
  endif
  r0 = r;
endfor

R3 = [0 ds0*cumsum(r.^(0:numel(R1)-2))]/L;

clear L r0 i r ds0 maxAR

% R3 is a much more agressive size distribution and that causes skewness as we
% transition from R1 to R3, therefore we make it milder at the cylinder corner.
w = 0.7;
R2 = (1-w)*R1+w*R3;

% mesh outer surface of sphere
Xan1 = zeros(numel(R1),data.N+1,Nz);
Yan1 = Xan1;  Zan1 = Xan1;

Xan1(1,:,:) = Xcyl(end,:,:);
Yan1(1,:,:) = Ycyl(end,:,:);
Zan1(1,:,:) = Zcyl(end,:,:);

Xan1(end,:,:) = Rsph * cos(phi)' * sin(theta);
Yan1(end,:,:) = Rsph * sin(phi)' * sin(theta);
Zan1(end,:,:) = Rsph * repmat(cos(theta),data.N+1,1);

% create inner layers by connecting the cylinder to the sphere
for i=1:Nz
  w = 1-((i-1)/(Nz-1))^1;
  w = w*R1' + (1-w)*R2';
  Xan1(:,:,i) = (1-w)*Xan1(1,:,i) + w*Xan1(end,:,i);
  Yan1(:,:,i) = (1-w)*Yan1(1,:,i) + w*Yan1(end,:,i);
  Zan1(:,:,i) = (1-w)*Zan1(1,:,i) + w*Zan1(end,:,i);
endfor

% Prolongate the anular domain around the cylinder corner, i.e. extend along the
% third axis, such that the base of the new section is the side of the cylinder.

Xan2 = zeros(numel(R1),data.N+1,data.M-BL+1);
Yan2 = Xan2;  Zan2 = Xan2;

Xan2(1,:,:) = fliplr(Xcyl(BL:end,:,end)');
Yan2(1,:,:) = fliplr(Ycyl(BL:end,:,end)');
Zan2(1,:,:) = fliplr(Zcyl(BL:end,:,end)');

% for this area new distributions of polar and azimuthal angle are required
% azimuthal angle taken from the cylinder end
phi = fliplr(atan2(Y(BL:end,:), X(BL:end,:))');

% polar angle set as scaled projection from a dummy focal point
focus = data.L*data.chord/tan(pi/4/data.projFactor);
r = fliplr(sqrt(X(BL:end,:).^2 + Y(BL:end,:).^2)');
theta = atan(r/focus)*data.projFactor;

Xan2(end,:,:) = Rsph * sin(theta).* cos(phi);
Yan2(end,:,:) = Rsph * sin(theta).* sin(phi);
Zan2(end,:,:) = Rsph * cos(theta);

clear r theta phi

N = data.M-BL+1;
for i=1:N
  w = ((N-i)/(N-1))^2;
  w = w*R2' + (1-w)*R3';
  Xan2(:,:,i) = (1-w)*Xan2(1,:,i) + w*Xan2(end,:,i);
  Yan2(:,:,i) = (1-w)*Yan2(1,:,i) + w*Yan2(end,:,i);
  Zan2(:,:,i) = (1-w)*Zan2(1,:,i) + w*Zan2(end,:,i);
endfor
clear i N w R1 R2

% merge anular shapes
Xan = cat(3,Xan1,Xan2(:,:,2:end));
Yan = cat(3,Yan1,Yan2(:,:,2:end));
Zan = cat(3,Zan1,Zan2(:,:,2:end));

clear Xan1 Xan2 Yan1 Yan2 Zan1 Zan2

% Blend edges (explicit smoothing)
radius = 0.4;

for i=1:4
  % 1st the corner of the anular shape
  for j=0:i
    for k=unique([Nz-j Nz+j])
      Xan(:,:,k) += 0.5*(Xan(:,:,k-1)+Xan(:,:,k+1)) - Xan(:,:,k);
      Yan(:,:,k) += 0.5*(Yan(:,:,k-1)+Yan(:,:,k+1)) - Yan(:,:,k);
      Zan(:,:,k) += 0.5*(Zan(:,:,k-1)+Zan(:,:,k+1)) - Zan(:,:,k);
    endfor
  endfor

  % update cylinder surface
  Xcyl(end,:,:) = Xan(1,:,1:Nz);
  Ycyl(end,:,:) = Yan(1,:,1:Nz);
  Zcyl(end,:,:) = Zan(1,:,1:Nz);

  Xcyl(BL:end,:,end) = fliplr(squeeze(Xan(1,:,Nz:end)))';
  Ycyl(BL:end,:,end) = fliplr(squeeze(Yan(1,:,Nz:end)))';
  Zcyl(BL:end,:,end) = fliplr(squeeze(Zan(1,:,Nz:end)))';

  % 2nd the top of the cylinder with the anular shape
  % j=0 (interface)
  U = 0.5*(Xcyl(end-1,:,:)+Xan(2,:,1:Nz)) - Xcyl(end,:,:);
  V = 0.5*(Ycyl(end-1,:,:)+Yan(2,:,1:Nz)) - Ycyl(end,:,:);
  W = 0.5*(Zcyl(end-1,:,:)+Zan(2,:,1:Nz)) - Zcyl(end,:,:);
  
  Xcyl(end,:,:) += U;  Xan(1,:,1:Nz) += U;
  Ycyl(end,:,:) += V;  Yan(1,:,1:Nz) += V;
  Zcyl(end,:,:) += W;  Zan(1,:,1:Nz) += W;
  
  for j=1:i
    % cylinder side
    Xcyl(end-j,:,:) += ...
      0.5*(Xcyl(end-j+1,:,:)+Xcyl(end-j-1,:,:)) - Xcyl(end-j,:,:);
    Ycyl(end-j,:,:) += ...
      0.5*(Ycyl(end-j+1,:,:)+Ycyl(end-j-1,:,:)) - Ycyl(end-j,:,:);
    Zcyl(end-j,:,:) += ...
      0.5*(Zcyl(end-j+1,:,:)+Zcyl(end-j-1,:,:)) - Zcyl(end-j,:,:);
    
    % annulus side
    Xan(1+j,:,1:Nz) += 0.5*(Xan(j,:,1:Nz)+Xan(2+j,:,1:Nz))-Xan(1+j,:,1:Nz);
    Yan(1+j,:,1:Nz) += 0.5*(Yan(j,:,1:Nz)+Yan(2+j,:,1:Nz))-Yan(1+j,:,1:Nz);
    Zan(1+j,:,1:Nz) += 0.5*(Zan(j,:,1:Nz)+Zan(2+j,:,1:Nz))-Zan(1+j,:,1:Nz);
  endfor

  % 3rd the end of the cylinder with the annulus, fixing the corner and the wall
  [~,~,N] = size(Xan(:,:,Nz:end));
  w = (0:N-1)/(N-1)/radius;
  w = 1 - max(1-w,0).^4.*(1+4*w) - max(1-fliplr(w),0).^4.*(1+4*fliplr(w));
  w = repmat(w,columns(Xan),1);
  w = reshape(w,[1 size(w)]);

  % j=0 (interface)
  U = w.*(0.5*(reshape(fliplr(Xcyl(BL:end,:,end-1)'),1,data.N+1,N)+...
      Xan(2,:,Nz:end)) - Xan(1,:,Nz:end));
  V = w.*(0.5*(reshape(fliplr(Ycyl(BL:end,:,end-1)'),1,data.N+1,N)+...
      Yan(2,:,Nz:end)) - Yan(1,:,Nz:end));
  W = w.*(0.5*(reshape(fliplr(Zcyl(BL:end,:,end-1)'),1,data.N+1,N)+...
      Zan(2,:,Nz:end)) - Zan(1,:,Nz:end));
  
  Xcyl(BL:end,:,end) += fliplr(squeeze(U))';  Xan(1,:,Nz:end) += U;
  Ycyl(BL:end,:,end) += fliplr(squeeze(V))';  Yan(1,:,Nz:end) += V;
  Zcyl(BL:end,:,end) += fliplr(squeeze(W))';  Zan(1,:,Nz:end) += W;
  
  for j=1:i
    % cylinder side
    Xcyl(BL:end,:,end-j) += squeeze(w)'.*(0.5*(Xcyl(BL:end,:,end-j+1) + ...
                            Xcyl(BL:end,:,end-j-1)) - Xcyl(BL:end,:,end-j));
    Ycyl(BL:end,:,end-j) += squeeze(w)'.*(0.5*(Ycyl(BL:end,:,end-j+1) + ...
                            Ycyl(BL:end,:,end-j-1)) - Ycyl(BL:end,:,end-j));
    Zcyl(BL:end,:,end-j) += squeeze(w)'.*(0.5*(Zcyl(BL:end,:,end-j+1) + ...
                            Zcyl(BL:end,:,end-j-1)) - Zcyl(BL:end,:,end-j));
    
    % annulus side
    Xan(1+j,:,Nz:end) += ...
      w.*(0.5*(Xan(j,:,Nz:end)+Xan(2+j,:,Nz:end))-Xan(1+j,:,Nz:end));
    Yan(1+j,:,Nz:end) += ...
      w.*(0.5*(Yan(j,:,Nz:end)+Yan(2+j,:,Nz:end))-Yan(1+j,:,Nz:end));
    Zan(1+j,:,Nz:end) += ...
      w.*(0.5*(Zan(j,:,Nz:end)+Zan(2+j,:,Nz:end))-Zan(1+j,:,Nz:end));
  endfor

  % 4th the corner, using 3-point averaging for it and 4-point on neighbours
  % corner
  Xan(1,:,Nz) = 1/3 * (Xan(1,:,Nz-1) + Xan(1,:,Nz+1) + Xan(2,:,Nz));
  Yan(1,:,Nz) = 1/3 * (Yan(1,:,Nz-1) + Yan(1,:,Nz+1) + Yan(2,:,Nz));
  Zan(1,:,Nz) = 1/3 * (Zan(1,:,Nz-1) + Zan(1,:,Nz+1) + Zan(2,:,Nz));
  
  for j=1:i
    % edge 1
    w = 3;
    Xan(1+j,:,Nz) = (Xan(j,:,Nz-1)+Xan(j,:,Nz+1)+Xan(j,:,Nz)+...
                     w*Xan(2+j,:,Nz))/(3+w);
    Yan(1+j,:,Nz) = (Yan(j,:,Nz-1)+Yan(j,:,Nz+1)+Yan(j,:,Nz)+...
                     w*Yan(2+j,:,Nz))/(3+w);
    Zan(1+j,:,Nz) = (Zan(j,:,Nz-1)+Zan(j,:,Nz+1)+Zan(j,:,Nz)+...
                     w*Zan(2+j,:,Nz))/(3+w);
    % edge 2
    w = 0.75;
    Xan(1,:,Nz+j) = (Xan(1,:,Nz+j-1)+Xan(2,:,Nz+j)+Xcyl(end-j,:,end-1)+...
                     w*Xan(1,:,Nz+j+1))/(3+w);
    Yan(1,:,Nz+j) = (Yan(1,:,Nz+j-1)+Yan(2,:,Nz+j)+Ycyl(end-j,:,end-1)+...
                     w*Yan(1,:,Nz+j+1))/(3+w);
    Zan(1,:,Nz+j) = (Zan(1,:,Nz+j-1)+Zan(2,:,Nz+j)+Zcyl(end-j,:,end-1)+...
                     w*Zan(1,:,Nz+j+1))/(3+w);
    % edge 3
    w = 1;
    Xan(1,:,Nz-j) = (Xan(1,:,Nz-j+1)+Xan(2,:,Nz-j)+Xcyl(end-1,:,end-j)+...
                     w*Xan(1,:,Nz-j-1))/(3+w);
    Yan(1,:,Nz-j) = (Yan(1,:,Nz-j+1)+Yan(2,:,Nz-j)+Ycyl(end-1,:,end-j)+...
                     w*Yan(1,:,Nz-j-1))/(3+w);
    Zan(1,:,Nz-j) = (Zan(1,:,Nz-j+1)+Zan(2,:,Nz-j)+Zcyl(end-1,:,end-j)+...
                     w*Zan(1,:,Nz-j-1))/(3+w);
  endfor

  % update cylinder surface
  Xcyl(end,:,:) = Xan(1,:,1:Nz);
  Ycyl(end,:,:) = Yan(1,:,1:Nz);
  Zcyl(end,:,:) = Zan(1,:,1:Nz);

  Xcyl(BL:end,:,end) = fliplr(squeeze(Xan(1,:,Nz:end)))';
  Ycyl(BL:end,:,end) = fliplr(squeeze(Yan(1,:,Nz:end)))';
  Zcyl(BL:end,:,end) = fliplr(squeeze(Zan(1,:,Nz:end)))';

endfor

%[~,~,N] = size(Xan);
%keep = ones(1,N);
%del = 2:2:Nz;
%keep(del) = 0;
%
%Xan = Xan(:,:,keep==1);
%Yan = Yan(:,:,keep==1);
%Zan = Zan(:,:,keep==1);
%
%Xcyl = Xcyl(1:end-1,:,:);
%Ycyl = Ycyl(1:end-1,:,:);
%Zcyl = Zcyl(1:end-1,:,:);

clear radius w i j k U V W N keep del
