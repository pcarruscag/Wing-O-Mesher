%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
%% Create a mesh for the solid domain %%
% Extrude the root section and apply the deformations applied to the fluid

printf('  Extruding base section...\n'); fflush(stdout);

z = surfSplits.zdist;

Xq0_s = repmat(section.xq0,1,1,Nz_s);
Yq0_s = repmat(section.yq0,1,1,Nz_s);
Zq0_s = repmat(z,size(section.xq0));

Xq1_s = repmat(section.xq1,1,1,Nz_s);
Yq1_s = repmat(section.yq1,1,1,Nz_s);
Zq1_s = repmat(z,size(section.xq1));

Xq2_s = repmat(section.xq2,1,1,Nz_s);
Yq2_s = repmat(section.yq2,1,1,Nz_s);
Zq2_s = repmat(z,size(section.xq2));

Xq3_s = repmat(section.xq3,1,1,Nz_s);
Yq3_s = repmat(section.yq3,1,1,Nz_s);
Zq3_s = repmat(z,size(section.xq3));

Xq4_s = repmat(section.xq4,1,1,Nz_s);
Yq4_s = repmat(section.yq4,1,1,Nz_s);
Zq4_s = repmat(z,size(section.xq4));

Xq5_s = repmat(section.xq5,1,1,Nz_s);
Yq5_s = repmat(section.yq5,1,1,Nz_s);
Zq5_s = repmat(z,size(section.xq5));

z = squeeze(z)';

Xt0_s = repmat(section.xt0,1,Nz_s);
Yt0_s = repmat(section.yt0,1,Nz_s);
Zt0_s = repmat(z,rows(section.xt0),1);

Xt1_s = repmat(section.xt1,1,Nz_s);
Yt1_s = repmat(section.yt1,1,Nz_s);
Zt1_s = repmat(z,rows(section.xt1),1);

Xt2_s = repmat(section.xt2,1,Nz_s);
Yt2_s = repmat(section.yt2,1,Nz_s);
Zt2_s = repmat(z,rows(section.xt2),1);

clear z section

% Smooth corner and apply the sweep and taper to solid domain
if data.smoothTip || data.sweep!=0 || data.taper!=1
  printf('  3D RBF mesh deformation:\n'); fflush(stdout);
endif

for i=1:2
  if i==1
    if(!data.smoothTip) continue; endif
    fun = defXYZ;
    rbf = smoother;
  else
    if(data.sweep==0 && data.taper==1) break; endif
    fun = defXY;
    rbf = sweeper;
  endif
  [Xq0_s,Yq0_s,Zq0_s] = defHexBlock(Xq0_s,Yq0_s,Zq0_s,fun,rbf,Nz_s);
  [Xq1_s,Yq1_s,Zq1_s] = defHexBlock(Xq1_s,Yq1_s,Zq1_s,fun,rbf,Nz_s);
  [Xq2_s,Yq2_s,Zq2_s] = defHexBlock(Xq2_s,Yq2_s,Zq2_s,fun,rbf,Nz_s);
  [Xq3_s,Yq3_s,Zq3_s] = defHexBlock(Xq3_s,Yq3_s,Zq3_s,fun,rbf,Nz_s);
  [Xq4_s,Yq4_s,Zq4_s] = defHexBlock(Xq4_s,Yq4_s,Zq4_s,fun,rbf,Nz_s);
  [Xq5_s,Yq5_s,Zq5_s] = defHexBlock(Xq5_s,Yq5_s,Zq5_s,fun,rbf,Nz_s);

  [Xt0_s,Yt0_s,Zt0_s] = defPrismBlock(Xt0_s,Yt0_s,Zt0_s,fun,rbf,Nz_s);
  [Xt1_s,Yt1_s,Zt1_s] = defPrismBlock(Xt1_s,Yt1_s,Zt1_s,fun,rbf,Nz_s);
  [Xt2_s,Yt2_s,Zt2_s] = defPrismBlock(Xt2_s,Yt2_s,Zt2_s,fun,rbf,Nz_s);
endfor
clear i fun rbf smoother defXYZ sweeper defXY

if(!EXPORT) return; endif

% Number nodes, count them and the elements
NELEM = 0;  NPOIN = 0;

% quad block 1
Ne = numel(Xq0_s(2:end,2:end,2:end));
Np = numel(Xq0_s);
Iq0_s = reshape(0:Np-1, size(Xq0_s));
NELEM += Ne;
NPOIN += Np;

% quad block 2
Ne = numel(Xq1_s(2:end,2:end,2:end));
Np = numel(Xq1_s(:,2:end,:));
Iq1_s = NPOIN*ones(size(Xq1_s));
Iq1_s(1:columns(Iq0_s),1,:) = Iq0_s(1,:,:);
Iq1_s(columns(Iq0_s):end-columns(Iq0_s)+1,1,:) = Iq0_s(:,end,:);
Iq1_s(end-columns(Iq0_s)+1:end,1,:) = fliplr(Iq0_s(end,:,:));
Iq1_s(:,2:end,:) += reshape(0:Np-1, size(Xq1_s(:,2:end,:)));
NELEM += Ne;
NPOIN += Np;

% triangle block 1
Ne = rows(tri0)*(Nz_s-1);
Np = numel(Xt0_s)-numel(Xq1_s(:,end,:));
It0_s = NPOIN*ones(size(Xt0_s));
It0_s(1:rows(Xq1_s),:) = Iq1_s(:,end,:);
It0_s(rows(Xq1_s)+1:end,:) += ...
  reshape(0:Np-1, rows(It0_s)-rows(Xq1_s), Nz_s);
NELEM += Ne;
NPOIN += Np;

% quad block 3
Ne = numel(Xq2_s(2:end,2:end,2:end));
Np = numel(Xq2_s(:,2:end,:));
Iq2_s = NPOIN*ones(size(Xq2_s));
Iq2_s(:,1,:) = It0_s(end-rows(Xq2_s)+1:end,:);
Iq2_s(:,2:end,:) += reshape(0:Np-1, size(Xq2_s(:,2:end,:)));
NELEM += Ne;
NPOIN += Np;

% triangle block 2
Ne = rows(tri1)*(Nz_s-1);
Np = numel(Xt1_s)-numel(Xq2_s(:,end,:));
It1_s = NPOIN*ones(size(Xt1_s));
It1_s(1:rows(Xq2_s),:) = Iq2_s(:,end,:);
It1_s(rows(Xq2_s)+1:end,:) += ...
  reshape(0:Np-1, rows(It1_s)-rows(Xq2_s), Nz_s);
NELEM += Ne;
NPOIN += Np;

% quad block 4
Ne = numel(Xq3_s(2:end,2:end,2:end));
Np = numel(Xq3_s(:,2:end,:));
Iq3_s = NPOIN*ones(size(Xq3_s));
Iq3_s(:,1,:) = It1_s(end-rows(Xq3_s)+1:end,:);
Iq3_s(:,2:end,:) += reshape(0:Np-1, size(Xq3_s(:,2:end,:)));
NELEM += Ne;
NPOIN += Np;

% triangle block 3
Ne = rows(tri2)*(Nz_s-1);
Np = numel(Xt2_s)-numel(Xq3_s(:,end,:));
It2_s = NPOIN*ones(size(Xt2_s));
It2_s(1:rows(Xq3_s),:) = Iq3_s(:,end,:);
It2_s(rows(Xq3_s)+1:end,:) += ...
  reshape(0:Np-1, rows(It2_s)-rows(Xq3_s), Nz_s);
NELEM += Ne;
NPOIN += Np;

% quad block 5
Ne = numel(Xq4_s(2:end,2:end,2:end));
Np = numel(Xq4_s(:,2:end,:));
Iq4_s = NPOIN*ones(size(Xq4_s));
Iq4_s(:,1,:) = It2_s(end-rows(Xq4_s)+1:end,:);
Iq4_s(:,2:end,:) += reshape(0:Np-1, size(Xq4_s(:,2:end,:)));
NELEM += Ne;
NPOIN += Np;

% quad block 6
Ne = numel(Xq5_s(2:end,2:end,2:end));
Np = numel(Xq5_s(2:end-1,2:end,:));
Iq5_s = NPOIN*ones(size(Xq5_s));
N = columns(Iq5_s);
Iq5_s(1,:,:) = flipud(squeeze(Iq4_s(1:N,end,:)));
Iq5_s(:,1,:) = squeeze(Iq4_s(N:end-N+1,end,:));
Iq5_s(end,:,:) = squeeze(Iq4_s(end-N+1:end,end,:));
Iq5_s(2:end-1,2:end,:) += reshape(0:Np-1, size(Xq5_s(2:end-1,2:end,:)));
NELEM += Ne;
NPOIN += Np;

printf("  Writing solid mesh...\n");
printf("   Num. vertices: %i\n",NPOIN);
printf("   Num. elements: %i\n",NELEM);
fflush(stdout);

clear N Ne Np

%%% Write mesh file %%%

fid = fopen("solid.su2","w");
fprintf(fid,"NDIME= 3\n");

% elements
fprintf(fid,"NELEM= %i\n",NELEM);

offset = writeHexBlock(fid,Iq0_s,0,0);
offset = writeHexBlock(fid,Iq1_s,offset,0);
offset = writeWedgeBlock(fid,tri0,It0_s,offset);
offset = writeHexBlock(fid,Iq2_s,offset,0);
offset = writeWedgeBlock(fid,tri1,It1_s,offset);
offset = writeHexBlock(fid,Iq3_s,offset,0);
offset = writeWedgeBlock(fid,tri2,It2_s,offset);
offset = writeHexBlock(fid,Iq4_s,offset,0);
offset = writeHexBlock(fid,Iq5_s,offset,0);

% points
fprintf(fid,"NPOIN= %i\n",NPOIN);

offset = write3Dgrid(fid,Xq0_s,Yq0_s,Zq0_s,[0 0 0 0 0 0],0);
offset = write3Dgrid(fid,Xq1_s,Yq1_s,Zq1_s,[0 0 1 0 0 0],offset);
mask = zeros(rows(Xt0_s),1);  mask(1:rows(Xq1_s)) = 1;
offset = writeWedgePts(fid,Xt0_s,Yt0_s,Zt0_s,mask,offset);
offset = write3Dgrid(fid,Xq2_s,Yq2_s,Zq2_s,[0 0 1 0 0 0],offset);
mask = zeros(rows(Xt1_s),1);  mask(1:rows(Xq2_s)) = 1;
offset = writeWedgePts(fid,Xt1_s,Yt1_s,Zt1_s,mask,offset);
offset = write3Dgrid(fid,Xq3_s,Yq3_s,Zq3_s,[0 0 1 0 0 0],offset);
mask = zeros(rows(Xt2_s),1);  mask(1:rows(Xq3_s)) = 1;
offset = writeWedgePts(fid,Xt2_s,Yt2_s,Zt2_s,mask,offset);
offset = write3Dgrid(fid,Xq4_s,Yq4_s,Zq4_s,[0 0 1 0 0 0],offset);
offset = write3Dgrid(fid,Xq5_s,Yq5_s,Zq5_s,[1 1 1 0 0 0],offset);

% boundaries
fprintf(fid,"NMARK= %i\n",2+5*surfSplits.number);

Ntip = NELEM/(Nz_s-1);
Nwing = (data.N-1)*(Nz_s-1);

fprintf(fid,"MARKER_TAG= clamped_s\n");
fprintf(fid,"MARKER_ELEMS= %i\n",Ntip);
writeIJboundary(fid,Iq0_s,0,1,true);
writeIJboundary(fid,Iq1_s,0,1,true);
writeTriPatch(fid,tri0,It0_s,1);
writeIJboundary(fid,Iq2_s,0,1,true);
writeTriPatch(fid,tri1,It1_s,1);
writeIJboundary(fid,Iq3_s,0,1,true);
writeTriPatch(fid,tri2,It2_s,1);
writeIJboundary(fid,Iq4_s,0,1,true);
writeIJboundary(fid,Iq5_s,0,1,true);

for i=1:surfSplits.number
  kRange = surfSplits.sections_s(i):surfSplits.sections_s(i+1);
  sk = numel(kRange)-1;

  % front part of suction side
  fprintf(fid,"MARKER_TAG= ss_f_%i_s\n",i);
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.nssf*sk);
  writeIKboundary(fid,Iq0_s,1,1:ceil(rows(Xq0_s)/2),kRange);
  writeJKboundary(fid,Iq1_s,0,rows(Iq1_s),[],kRange);
  mask = mask0;  mask(2:rows(Iq1_s)-1) = 0;
  writeWedgeEnds(fid,It0_s,mask,kRange,false,true);
  writeJKboundary(fid,Iq2_s,0,rows(Iq2_s),[],kRange);
  mask = mask1;  mask(2:rows(Iq2_s)-1) = 0;
  writeWedgeEnds(fid,It1_s,mask,kRange,false,true);

  % back part of suction side
  fprintf(fid,"MARKER_TAG= ss_b_%i_s\n",i);
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.nssb*sk);
  writeJKboundary(fid,Iq3_s,0,rows(Iq3_s),[],kRange);

  % trailing edge
  fprintf(fid,"MARKER_TAG= te_%i_s\n",i);
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.nte*sk);
  mask = mask2;  mask(2:rows(Iq3_s)-1) = 0;
  writeWedgeEnds(fid,It2_s,mask,kRange,true,false,true);
  writeWedgeEnds(fid,It2_s,mask,kRange,false,true);
  writeJKboundary(fid,Iq4_s,0,1,[],kRange,true);
  writeJKboundary(fid,Iq4_s,0,rows(Iq4_s),[],kRange);
  writeIKboundary(fid,Iq5_s,columns(Iq5_s),[],kRange,true);

  % back part of pressure side
  fprintf(fid,"MARKER_TAG= ps_b_%i_s\n",i);
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.npsb*sk);
  writeJKboundary(fid,Iq3_s,0,1,[],kRange,true);

  % front part of pressure side
  fprintf(fid,"MARKER_TAG= ps_f_%i_s\n",i);
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.npsf*sk);
  writeIKboundary(fid,Iq0_s,1,ceil(rows(Xq0_s)/2):rows(Xq0_s),kRange);
  writeJKboundary(fid,Iq1_s,0,1,[],kRange,true);
  mask = mask0;  mask(2:rows(Iq1_s)-1) = 0;
  writeWedgeEnds(fid,It0_s,mask,kRange,true,false,true);
  writeJKboundary(fid,Iq2_s,0,1,[],kRange,true);
  mask = mask1;  mask(2:rows(Iq2_s)-1) = 0;
  writeWedgeEnds(fid,It1_s,mask,kRange,true,false,true);
endfor

fprintf(fid,"MARKER_TAG= tip_s\n");
fprintf(fid,"MARKER_ELEMS= %i\n",Ntip);
writeIJboundary(fid,Iq0_s,0,Nz_s);
writeIJboundary(fid,Iq1_s,0,Nz_s);
writeTriPatch(fid,tri0,It0_s,Nz_s,true);
writeIJboundary(fid,Iq2_s,0,Nz_s);
writeTriPatch(fid,tri1,It1_s,Nz_s,true);
writeIJboundary(fid,Iq3_s,0,Nz_s);
writeTriPatch(fid,tri2,It2_s,Nz_s,true);
writeIJboundary(fid,Iq4_s,0,Nz_s);
writeIJboundary(fid,Iq5_s,0,Nz_s);

fclose(fid);

clear fid offset ans mask Ntip Nwing NPOIN NELEM Nz_s i kRange sk
