%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
%%% Assign node numbers and count nodes and elements %%%
NELEM = 0;
NPOIN = 0;

% cylinder
[M,N,P] = size(Xcyl);
Icyl = zeros(M,N,P);
Np = M*(N-2)*P;
Ne = (M-1)*(N-2)*(P-1);
Icyl(:,2:end-1,:) = reshape(0:Np-1, M, N-2, P);
Icyl(:,end,:) = Icyl(:,2,:);
Icyl(:,1,:) = Icyl(:,end-1,:);
NELEM += Ne;
NPOIN += Np;

% annulus
[M,N,P] = size(Xan);
Ian = ones(M,N,P)*NPOIN;
Np = (M-1)*(N-2)*P;
Ne = (M-1)*(N-2)*(P-1);
Ian(1,:,1:Nz) = Icyl(end,:,:);
Ian(1,:,Nz:end) = fliplr(Icyl(BL:end,:,end)');
Ian(2:end,2:end-1,:) += reshape(0:Np-1, M-1, N-2, P);
Ian(:,end,:) = Ian(:,2,:);
Ian(:,1,:) = Ian(:,end-1,:);
NELEM += Ne;
NPOIN += Np;

% tip extrusion
% quad block 1
[M,N,P] = size(Xq0);
Iq0 = ones(M,N,P)*NPOIN;
Np = M*(N-1)*P;
Ne = (M-1)*(N-1)*(P-1);
mid = ceil(M/2);
Iq0(1:mid,1,1:BL) = Icyl(1:BL,end-mid+1:end,end)';
Iq0(mid:end,1,1:BL) = Icyl(1:BL,2:mid+1,end)';
Iq0(1:mid,1,BL:end) = Ian(:,end-mid+1:end,end)';
Iq0(mid:end,1,BL:end) = Ian(:,2:mid+1,end)';
Iq0(:,2:end,:) += reshape(0:Np-1, M, N-1, P);
NELEM += Ne;
NPOIN += Np;

% quad block 2
[M,N,P] = size(Xq1);
Iq1 = ones(M,N,P)*NPOIN;
Np = (M-2)*(N-1)*P;
Ne = (M-1)*(N-1)*(P-1);
Iq1(1:columns(Iq0),1,:) = Iq0(1,:,:);
Iq1(columns(Iq0):end-columns(Iq0)+1,1,:) = Iq0(:,end,:);
Iq1(end-columns(Iq0)+1:end,1,:) = fliplr(Iq0(end,:,:));
Iq1(1,:,1:BL) = flipud(Icyl(1:BL,cut2+1:end-mid+1,end)');
Iq1(1,:,BL:end) = flipud(Ian(:,cut2+1:end-mid+1,end)');
Iq1(end,:,1:BL) = Icyl(1:BL,mid+1:cut1+1,end)';
Iq1(end,:,BL:end) = Ian(:,mid+1:cut1+1,end)';
Iq1(2:end-1,2:end,:) += reshape(0:Np-1, M-2, N-1, P);
NELEM += Ne;
NPOIN += Np;

% triangle block 1
[N,P] = size(Xt0);
It0 = ones(N,P)*NPOIN;
mask0 = ones(M,1);
for i=1:cut3-cut1
  mask0 = [mask0; cut2+1-i; zeros(max(M-2*i-2,nLE-1),1); cut1+1+i];
endfor
Np = sum(mask0==0)*P;
Ne = rows(tri0)*(P-1);
k = 0;
for i=1:P
  It0(1:M,i) = Iq1(:,end,i);
  
  for j=M+1:N
    if mask0(j)==0
      It0(j,i) += k++;
    else
      if i<BL
        It0(j,i) = Icyl(i,mask0(j),end);
      else
        It0(j,i) = Ian(i-BL+1,mask0(j),end);
      endif
    endif
  endfor
endfor
NELEM += Ne;
NPOIN += Np;

% quad block 3
[M,N,P] = size(Xq2);
Iq2 = ones(M,N,P)*NPOIN;
Np = (M-2)*(N-1)*P;
Ne = (M-1)*(N-1)*(P-1);
for i=1:P
  Iq2(:,1,i) = It0(end-M+1:end,i);
endfor
Iq2(1,:,1:BL) = flipud(Icyl(1:BL,cut6+1:cut4+1,end)');
Iq2(1,:,BL:end) = flipud(Ian(:,cut6+1:cut4+1,end)');
Iq2(end,:,1:BL) = Icyl(1:BL,cut3+1:cut5+1,end)';
Iq2(end,:,BL:end) = Ian(:,cut3+1:cut5+1,end)';
Iq2(2:end-1,2:end,:) += reshape(0:Np-1, M-2, N-1, P);
NELEM += Ne;
NPOIN += Np;

% triangle block 2
[N,P] = size(Xt1);
It1 = ones(N,P)*NPOIN;
mask1 = ones(M,1);
for i=1:cut7-cut5
  mask1 = [mask1; cut6+1-i; zeros(max(M-i-2,nTAR-1),1); cut5+1+i];
endfor
Np = sum(mask1==0)*P;
Ne = rows(tri1)*(P-1);
k = 0;
for i=1:P
  It1(1:M,i) = Iq2(:,end,i);
  
  for j=M+1:N
    if mask1(j)==0
      It1(j,i) += k++;
    else
      if i<BL
        It1(j,i) = Icyl(i,mask1(j),end);
      else
        It1(j,i) = Ian(i-BL+1,mask1(j),end);
      endif
    endif
  endfor
endfor
NELEM += Ne;
NPOIN += Np;

% quad block 4
[M,N,P] = size(Xq3);
Iq3 = ones(M,N,P)*NPOIN;
Np = (M-2)*(N-1)*P;
Ne = (M-1)*(N-1)*(P-1);
for i=1:P
  Iq3(:,1,i) = It1(end-M+1:end,i);
endfor
Iq3(1,:,1:BL) = flipud(Icyl(1:BL,cut10+1:cut8+1,end)');
Iq3(1,:,BL:end) = flipud(Ian(:,cut10+1:cut8+1,end)');
Iq3(end,:,1:BL) = Icyl(1:BL,cut7+1:cut9+1,end)';
Iq3(end,:,BL:end) = Ian(:,cut7+1:cut9+1,end)';
Iq3(2:end-1,2:end,:) += reshape(0:Np-1, M-2, N-1, P);
NELEM += Ne;
NPOIN += Np;

% triangle block 3
[N,P] = size(Xt2);
It2 = ones(N,P)*NPOIN;
mask2 = [];
for i=0:cut11-cut9
  mask2 = [cut12+1+i; zeros(max(rows(Xq4)-2*i-2,nTAR-1),1); cut11+1-i; mask2];
endfor
mask2(1:nTAR+1) = 1;
Np = sum(mask2==0)*P;
Ne = rows(tri2)*(P-1);
k = 0;
for i=1:P
  It2(1:M,i) = Iq3(:,end,i);
  
  for j=M+1:N
    if mask2(j)==0
      It2(j,i) += k++;
    else
      if i<BL
        It2(j,i) = Icyl(i,mask2(j),end);
      else
        It2(j,i) = Ian(i-BL+1,mask2(j),end);
      endif
    endif
  endfor
endfor
NELEM += Ne;
NPOIN += Np;

% quad block 5
[M,N,P] = size(Xq4);
Iq4 = ones(M,N,P)*NPOIN;
Np = (M-2)*(N-1)*P;
Ne = (M-1)*(N-1)*(P-1);
cut13 = cut11+N-1;
cut14 = cut12-N+1;
for i=1:P
  Iq4(:,1,i) = It2(end-M+1:end,i);
endfor
Iq4(1,:,1:BL) = flipud(Icyl(1:BL,cut14+1:cut12+1,end)');
Iq4(1,:,BL:end) = flipud(Ian(:,cut14+1:cut12+1,end)');
Iq4(end,:,1:BL) = Icyl(1:BL,cut11+1:cut13+1,end)';
Iq4(end,:,BL:end) = Ian(:,cut11+1:cut13+1,end)';
Iq4(2:end-1,2:end,:) += reshape(0:Np-1, M-2, N-1, P);
NELEM += Ne;
NPOIN += Np;

% quad block 6
[M,N,P] = size(Xq5);
Iq5 = ones(M,N,P)*NPOIN;
Np = (M-2)*(N-2)*P;
Ne = (M-1)*(N-1)*(P-1);
Iq5(1,:,:) = flipud(squeeze(Iq4(1:N,end,:)));
Iq5(:,1,:) = squeeze(Iq4(N:end-N+1,end,:));
Iq5(end,:,:) = squeeze(Iq4(end-N+1:end,end,:));
Iq5(:,end,1:BL) = flipud(Icyl(1:BL,cut13+1:cut14+1,end)');
Iq5(:,end,BL:end) = flipud(Ian(:,cut13+1:cut14+1,end)');
Iq5(2:end-1,2:end-1,:) += reshape(0:Np-1, M-2, N-2, P);
NELEM += Ne;
NPOIN += Np;

printf("  Writing fluid mesh...\n");
printf("   Num. vertices: %i\n",NPOIN);
printf("   Num. elements: %i\n",NELEM);
fflush(stdout);


%%% Write mesh file %%%

fid = fopen("mesh.su2","w");
fprintf(fid,"NDIME= 3\n");

%% elements %%
fprintf(fid,"NELEM= %i\n",NELEM);

% cylinder
offset = writeHexBlock(fid,Icyl,0,1);

% annulus
offset = writeHexBlock(fid,Ian,offset,1);

% tip extrusion
% quad block 1
offset = writeHexBlock(fid,Iq0,offset,0);

% quad block 2
offset = writeHexBlock(fid,Iq1,offset,0);

% triangle block 1
offset = writeWedgeBlock(fid,tri0,It0,offset);

% quad block 3
offset = writeHexBlock(fid,Iq2,offset,0);

% triangle block 2
offset = writeWedgeBlock(fid,tri1,It1,offset);

% quad block 4
offset = writeHexBlock(fid,Iq3,offset,0);

% triangle block 3
offset = writeWedgeBlock(fid,tri2,It2,offset);

% quad block 5
offset = writeHexBlock(fid,Iq4,offset,0);

% quad block 6
offset = writeHexBlock(fid,Iq5,offset,0);

%% points %%
fprintf(fid,"NPOIN= %i\n",NPOIN);

% cylinder
offset = write3Dgrid(fid,Xcyl,Ycyl,Zcyl,[0 0 1 1 0 0],0);

% annulus
offset = write3Dgrid(fid,Xan,Yan,Zan,[1 0 1 1 0 0],offset);

% tip extrusion
% quad block 1
offset = write3Dgrid(fid,Xq0,Yq0,Zq0,[0 0 1 0 0 0],offset);

% quad block 2
offset = write3Dgrid(fid,Xq1,Yq1,Zq1,[1 1 1 0 0 0],offset);

% triangle block 1
offset = writeWedgePts(fid,Xt0,Yt0,Zt0,mask0,offset);

% quad block 3
offset = write3Dgrid(fid,Xq2,Yq2,Zq2,[1 1 1 0 0 0],offset);

% triangle block 2
offset = writeWedgePts(fid,Xt1,Yt1,Zt1,mask1,offset);

% quad block 4
offset = write3Dgrid(fid,Xq3,Yq3,Zq3,[1 1 1 0 0 0],offset);

% triangle block 3
offset = writeWedgePts(fid,Xt2,Yt2,Zt2,mask2,offset);

% quad block 5
offset = write3Dgrid(fid,Xq4,Yq4,Zq4,[1 1 1 0 0 0],offset);

% quad block 6
offset = write3Dgrid(fid,Xq5,Yq5,Zq5,[1 1 1 1 0 0],offset);

%% boundaries %%
fprintf(fid,"NMARK= %i\n",3+5*surfSplits.number);

Ntip = rows(tri0) + rows(tri1) + rows(tri2) + ...
       numel(Iq0(2:end,2:end,1)) + numel(Iq1(2:end,2:end,1)) + ...
       numel(Iq2(2:end,2:end,1)) + numel(Iq3(2:end,2:end,1)) + ...
       numel(Iq4(2:end,2:end,1)) + numel(Iq5(2:end,2:end,1));

fprintf(fid,"MARKER_TAG= farfield\n");
[M,N,P] = size(Ian);
fprintf(fid,"MARKER_ELEMS= %i\n",(N-2)*(P-1)+Ntip);
writeJKboundary(fid,Ian,1,M);
M = columns(It0);
writeIJboundary(fid,Iq0,0,M);
writeIJboundary(fid,Iq1,0,M);
writeTriPatch(fid,tri0,It0,M,true);
writeIJboundary(fid,Iq2,0,M);
writeTriPatch(fid,tri1,It1,M,true);
writeIJboundary(fid,Iq3,0,M);
writeTriPatch(fid,tri2,It2,M,true);
writeIJboundary(fid,Iq4,0,M);
writeIJboundary(fid,Iq5,0,M);

fprintf(fid,"MARKER_TAG= symmetry\n");
[M1,N,~] = size(Icyl);
[M2,~,~] = size(Ian);
fprintf(fid,"MARKER_ELEMS= %i\n",(N-2)*(M1+M2-2));
writeIJboundary(fid,Icyl,1,1,true);
writeIJboundary(fid,Ian,1,1,true);

for i=1:surfSplits.number
  kRange = surfSplits.sections_f(i):surfSplits.sections_f(i+1);
  sk = numel(kRange)-1;

  % front part of suction side
  fprintf(fid,"MARKER_TAG= ss_f_%i\n",i);
  jRange = 1:cut7;
  surfSplits.nssf = numel(jRange)-1;
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.nssf*sk);
  writeJKboundary(fid,Icyl,1,1,jRange,kRange,true);

  % back part of suction side
  fprintf(fid,"MARKER_TAG= ss_b_%i\n",i);
  jRange = cut7:cut9;
  surfSplits.nssb = numel(jRange)-1;
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.nssb*sk);
  writeJKboundary(fid,Icyl,1,1,jRange,kRange,true);

  % trailing edge
  fprintf(fid,"MARKER_TAG= te_%i\n",i);
  jRange = cut9:cut10;
  surfSplits.nte = numel(jRange)-1;
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.nte*sk);
  writeJKboundary(fid,Icyl,1,1,jRange,kRange,true);

  % back part of pressure side
  fprintf(fid,"MARKER_TAG= ps_b_%i\n",i);
  jRange = cut10:cut8;
  surfSplits.npsb = numel(jRange)-1;
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.npsb*sk);
  writeJKboundary(fid,Icyl,1,1,jRange,kRange,true);

  % front part of pressure side
  fprintf(fid,"MARKER_TAG= ps_f_%i\n",i);
  jRange = cut8:columns(Icyl)-1;
  surfSplits.npsf = numel(jRange)-1;
  fprintf(fid,"MARKER_ELEMS= %i\n",surfSplits.npsf*sk);
  writeJKboundary(fid,Icyl,1,1,jRange,kRange,true);
endfor

fprintf(fid,"MARKER_TAG= tip\n");
fprintf(fid,"MARKER_ELEMS= %i\n",Ntip);
writeIJboundary(fid,Iq0,0,1,true);
writeIJboundary(fid,Iq1,0,1,true);
writeTriPatch(fid,tri0,It0,1);
writeIJboundary(fid,Iq2,0,1,true);
writeTriPatch(fid,tri1,It1,1);
writeIJboundary(fid,Iq3,0,1,true);
writeTriPatch(fid,tri2,It2,1);
writeIJboundary(fid,Iq4,0,1,true);
writeIJboundary(fid,Iq5,0,1,true);

fclose(fid);

clear M M1 M2 N P NELEM NPOIN Ne Np Ntip offset
clear i j k cut* nLE nTAR fid mid BL Nz ans sk jRange kRange
