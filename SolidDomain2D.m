%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
%% Create a mesh for the solid domain %%
assert(data.wake > 0, "Incompatible wake settings for solid mesh.");
if data.wake == 1
  nHalfThick = 2^data.refineFactor;
elseif data.wake == 2
  nHalfThick = 2*(data.refineFactor+1);
endif

% location along surface to start the domain, cannot be LE!
start = data.splitSurfaceAt;

[~,N] = size(X);
midPt = N/2+1;
Ms = nHalfThick*2+1;
Ns = midPt-nHalfThick-start;

% get surface nodes from fluid mesh
Xs = zeros(Ms,Ns); Ys = Xs;
Xs(1,:) = X(1,(1+start):(midPt-nHalfThick));
Ys(1,:) = Y(1,(1+start):(midPt-nHalfThick));
Xs(end,:) = fliplr(X(1,(midPt+nHalfThick):(N-start+1)));
Ys(end,:) = fliplr(Y(1,(midPt+nHalfThick):(N-start+1)));
clear N

% linear interpolation to create nodes inside the domain
w = (0:2*nHalfThick)'/nHalfThick/2;
Xs = (1-w)*Xs(1,:)+w*Xs(end,:);
Ys = (1-w)*Ys(1,:)+w*Ys(end,:);
clear w

% use fluid coordinates for the trailing edge
Xs(2:end-1,end) = X(1,(-nHalfThick+1:nHalfThick-1)+midPt);
Ys(2:end-1,end) = Y(1,(-nHalfThick+1:nHalfThick-1)+midPt);


% create more structural nodes
if data.factorM != 1 || data.factorN != 1
  newNs = ceil((Ns-1)*data.factorN);
  % "cofactorize" in the N direction
  while gcd(Ns,newNs++) > 1
  endwhile
  newMs = ceil((Ms-1)*data.factorM/2)*2+1;
  Xnew = zeros(newMs,newNs); Ynew = Xnew;
  ref = (0:Ns-1)/(Ns-1);
  new = (0:newNs-1)/(newNs-1);
  Xnew(1,:) = interp1(ref,Xs(1,:),new);
  Ynew(1,:) = interp1(ref,Ys(1,:),new);
  Xnew(end,:) = interp1(ref,Xs(end,:),new);
  Ynew(end,:) = interp1(ref,Ys(end,:),new);

  w = (0:newMs-1)'/(newMs-1);
  Xnew = (1-w)*Xnew(1,:)+w*Xnew(end,:);
  Ynew = (1-w)*Ynew(1,:)+w*Ynew(end,:);

  ref = (0:Ms-1)/(Ms-1);
  Xnew(:,end) = interp1(ref,Xs(:,end),w);
  Ynew(:,end) = interp1(ref,Ys(:,end),w);

  Xs = Xnew;  Ys = Ynew;
  Ns = newNs; Ms = newMs;
  nHalfThick = (Ms-1)/2;
  clear w ref Xnew Ynew newNs newMs new
endif


% the spacing of the nodes at the TE is far from linear, so we smooth near it
blendDist = 5;
idx = (Ns-blendDist):Ns;
s1 = (0:2*nHalfThick)'/nHalfThick/2;
s2 = [0; cumsum(sqrt(diff(Xs(:,end)).^2+diff(Ys(:,end)).^2))]; s2 /= s2(end);
w = (0:blendDist)/blendDist; w = s1*(1-w)+s2*w;
Xs(:,idx) = (1-w).*Xs(1,idx)+w.*Xs(end,idx);
Ys(:,idx) = (1-w).*Ys(1,idx)+w.*Ys(end,idx);
clear blendDist idx s1 s2 w


% refine
for k=1:data.refineFactorS
  % new sizes
  Ms = 2*Ms-1;  Ns = 2*Ns-1;
  nHalfThick = (Ms-1)/2;
  % interpolate
  isrc = 1:2:Ms;  jsrc = 1:2:Ns;
  [XI,YI] = meshgrid(1:Ns,1:Ms);
  Xs = interp2(jsrc,isrc,Xs,XI,YI,'linear');
  Ys = interp2(jsrc,isrc,Ys,XI,YI,'linear');
endfor
clear k isrc jsrc XI YI


% close the leading edge by extruding inward leaving a hole in the structure
if data.leLayers
  Mle = data.leLayers+1;
  Nle = 2*start-1;
  [~,N] = size(X);

  % outer nodes copied from fluid side
  Xle = zeros(Mle,Nle);  Yle = zeros(Mle,Nle);
  Xle(1,:) = [X(1,(N-start+1):end) X(1,3:(1+start))];
  Yle(1,:) = [Y(1,(N-start+1):end) Y(1,3:(1+start))];
  
  % the ends of this segment are set equal to the generated structural mesh
  Xle(:,end) = Xs(1:Mle,1);  Xle(:,1) = flipud(Xs((end-Mle+1):end,1));
  Yle(:,end) = Ys(1:Mle,1);  Yle(:,1) = flipud(Ys((end-Mle+1):end,1));

  % inner nodes by extrusion, at the LE thickness is reduced (self-intersect)
  tscale = data.leLayerThin(1)-...
           data.leLayerThin(2)*cos((-(start-2):(start-2))/(start-2)*pi/2);

  for i=2:Mle
    % compute normals
    Nx = Yle(i-1,3:end)-Yle(i-1,1:end-2);
    Ny = Xle(i-1,1:end-2)-Xle(i-1,3:end);
    
    % thickness
    t = 0.5*(sqrt((Xle(i,1)-Xle(i-1,1))^2+(Yle(i,1)-Yle(i-1,1))^2)+...
        sqrt((Xle(i,end)-Xle(i-1,end))^2+(Yle(i,end)-Yle(i-1,end))^2));
    t = t./sqrt(Nx.^2+Ny.^2);
    
    Xle(i,2:end-1) = Xle(i-1,2:end-1)+t.*tscale.*Nx;
    Yle(i,2:end-1) = Yle(i-1,2:end-1)+t.*tscale.*Ny;
  endfor
  clear i N t tscale Nx Ny
else
  Nle = 1;
endif

% plot
%clf
hold on
plot(Xs',Ys',"b")
plot(Xs,Ys,"b")
if data.leLayers
  plot(Xle',Yle',"b")
  plot(Xle,Yle,"b")
endif
axis equal
hold off

if(!EXPORT) return; endif

%% SAVE %%
fid = fopen("solid.su2","w");
fprintf(fid,"NDIME= 2\n");

% global vertex and element numbering
Nv = Ms*Ns;
Vidx = reshape((0:Nv-1)',Ns,Ms)';
Ne = (Ms-1)*(Ns-1);
Eidx = reshape((0:Ne-1)',Ns-1,Ms-1)';

if data.leLayers
  VidxLe = zeros(Mle,Nle);
  VidxLe(:,end) = Vidx(1:Mle,1);  VidxLe(:,1) = flipud(Vidx((end-Mle+1):end,1));
  for i=1:Mle
    for j=2:Nle-1
      VidxLe(i,j) = Nv++;
    endfor
  endfor

  EidxLe = zeros(Mle-1,Nle-1);
  for i=1:Mle-1
    for j=1:Nle-1
      EidxLe(i,j) = Ne++;
    endfor
  endfor
endif

printf("  Writing solid mesh...\n");
printf("   Num. vertices: %i\n",Nv);
printf("   Num. elements: %i\n",Ne);
fflush(stdout);

% element definition
fprintf(fid,"NELEM= %i\n",Ne);

function writeElems(fid,Eidx,Vidx)
  [M,N] = size(Eidx);
  for i=1:M
    for j=1:N
      % the element index
      eli = Eidx(i,j);
      % the vertices that define it
      v0 = Vidx( i , j );
      v1 = Vidx( i ,j+1);
      v2 = Vidx(i+1,j+1);
      v3 = Vidx(i+1, j );
      fprintf(fid,"9  %i  %i  %i  %i  %i\n",v0,v1,v2,v3,eli);
    endfor
  endfor
endfunction

writeElems(fid,Eidx,Vidx);

if data.leLayers
  writeElems(fid,EidxLe,VidxLe);
endif

% point coordinates
fprintf(fid,"NPOIN= %i\n",Nv);

for i=1:Ms
  for j=1:Ns
    fprintf(fid,"%e  %e  %i\n",Xs(i,j),Ys(i,j),Vidx(i,j));
  endfor
endfor

if data.leLayers
  for i=1:Mle
    for j=2:Nle-1
      fprintf(fid,"%e  %e  %i\n",Xle(i,j),Yle(i,j),VidxLe(i,j));
    endfor
  endfor
endif

% boundary definition
fprintf(fid,"NMARK= %i\n",3+(data.leLayers>0));

fprintf(fid,"MARKER_TAG= clamped\n");
fprintf(fid,"MARKER_ELEMS= %i\n",Ms-1-2*data.leLayers+Nle-1);

for i=(1+data.leLayers):(Ms-1-data.leLayers)
  fprintf(fid,"3  %i  %i\n",Vidx(i,1),Vidx(i+1,1));
endfor
if data.leLayers
  for i=1:Nle-1
    fprintf(fid,"3  %i  %i\n",VidxLe(end,i),VidxLe(end,i+1));
  endfor
endif

fprintf(fid,"MARKER_TAG= suction_side_s\n");
fprintf(fid,"MARKER_ELEMS= %i\n",Ns-1+nHalfThick);

for j=1:Ns-1
  fprintf(fid,"3  %i  %i\n",Vidx(1,j),Vidx(1,j+1));
endfor
for i=1:nHalfThick
  fprintf(fid,"3  %i  %i\n",Vidx(i,end),Vidx(i+1,end));
endfor

fprintf(fid,"MARKER_TAG= pressure_side_s\n");
fprintf(fid,"MARKER_ELEMS= %i\n",Ns-1+nHalfThick);

for i=1:nHalfThick
  fprintf(fid,"3  %i  %i\n",Vidx(nHalfThick+i,end),Vidx(nHalfThick+i+1,end));
endfor
for j=1:Ns-1
  fprintf(fid,"3  %i  %i\n",Vidx(end,Ns-j+1),Vidx(end,Ns-j));
endfor

if data.leLayers
  fprintf(fid,"MARKER_TAG= leading_edge_s\n");
  fprintf(fid,"MARKER_ELEMS= %i\n",Nle-1);

  for i=1:Nle-1
    fprintf(fid,"3  %i  %i\n",VidxLe(1,i),VidxLe(1,i+1));
  endfor
endif

fclose(fid);

clear Ne Nv ans fid i j nHalfThick Vidx Eidx start midPt Ms Ns
clear VidxLe EidxLe Mle Nle
