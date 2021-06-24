%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
1+1; clear ans % little hack to make this a "header" file

%% Write the hex elements of a 3D grid
function nel = writeHexBlock(fid,idx,offset,annulus)
  [M,N,P] = size(idx);
  v0 = idx(1:M-1, 1+annulus:N-1, 1:P-1);  v0 = reshape(v0,1,numel(v0));
  v1 = idx(1:M-1,  2+annulus:N,  1:P-1);  v1 = reshape(v1,1,numel(v1));
  v2 = idx( 2:M,   2+annulus:N,  1:P-1);  v2 = reshape(v2,1,numel(v2));
  v3 = idx( 2:M,  1+annulus:N-1, 1:P-1);  v3 = reshape(v3,1,numel(v3));
  v4 = idx(1:M-1, 1+annulus:N-1,  2:P );  v4 = reshape(v4,1,numel(v4));
  v5 = idx(1:M-1,  2+annulus:N,   2:P );  v5 = reshape(v5,1,numel(v5));
  v6 = idx( 2:M,   2+annulus:N,   2:P );  v6 = reshape(v6,1,numel(v6));
  v7 = idx( 2:M,  1+annulus:N-1,  2:P );  v7 = reshape(v7,1,numel(v7));
  eli = (0:((M-1)*(N-1-annulus)*(P-1)-1))+offset;
  fprintf(fid,"12  %i  %i  %i  %i  %i  %i  %i  %i  %i\n",...
                  [v0; v3; v2; v1; v4; v7; v6; v5; eli]);
  nel = eli(end)+1;
endfunction


%% Write wedge (triangular prism) elements
function nel = writeWedgeBlock(fid,tri,idx,offset)
  M = rows(tri);
  [N,P] = size(idx);

  for i=1:P-1
    v0 = idx(tri(:,1), i )';  v1 = idx(tri(:,2), i )';  v2 = idx(tri(:,3), i )';
    v3 = idx(tri(:,1),i+1)';  v4 = idx(tri(:,2),i+1)';  v5 = idx(tri(:,3),i+1)';
    eli = (0:M-1)+offset;
    fprintf(fid,"13  %i  %i  %i  %i  %i  %i  %i\n",...
                    [v0; v1; v2; v3; v4; v5; eli]);
    offset = eli(end)+1;
  endfor

  nel = offset;
endfunction


%% Write the point coordinates of a 3D hex grid
function npt = write3Dgrid(fid,X,Y,Z,dominated,offset)
  [M,N,P] = size(X);
  i = 1+dominated(1):M-dominated(2);
  j = 1+dominated(3):N-dominated(4);
  k = 1+dominated(5):P-dominated(6);
  npt = numel(i)*numel(j)*numel(k);
  x = reshape(X(i,j,k),1,npt);
  y = reshape(Y(i,j,k),1,npt);
  z = reshape(Z(i,j,k),1,npt);
  idx = (0:npt-1)+offset;
  fprintf(fid,"%e  %e  %e  %i\n",[x; z; y; idx]);
  npt = idx(end)+1;
endfunction


%% Write the point coordinates of a wedge block
function npt = writeWedgePts(fid,X,Y,Z,dominated,offset)
  [N,P] = size(X);
  keep = repmat(dominated==0,P,1);
  x = reshape(X,1,N*P)(keep);
  y = reshape(Y,1,N*P)(keep);
  z = reshape(Z,1,N*P)(keep);
  npt = sum(keep);
  idx = (0:npt-1)+offset;
  fprintf(fid,"%e  %e  %e  %i\n",[x; z; y; idx]);
  npt = idx(end)+1;
endfunction


%% Write IJ boundary of a 3D grid (quad elements)
function writeIJboundary(fid,idx,annulus,k,flip=false)
  M = numel(idx(2:end,2+annulus:end,1));
  v0 = reshape(idx(1:end-1, 1+annulus:end-1, k),1,M);
  v1 = reshape(idx(1:end-1,  2+annulus:end,  k),1,M);
  v2 = reshape(idx( 2:end,   2+annulus:end,  k),1,M);
  v3 = reshape(idx( 2:end,  1+annulus:end-1, k),1,M);
  if flip
    fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v3; v2; v1]);
  else
    fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v1; v2; v3]);
  endif
endfunction


%% Write JK boundary of a 3D grid (quad elements)
function writeJKboundary(fid,idx,annulus,i,jRange=[],kRange=[],flip=false)
  [~,N,P] = size(idx);
  if isempty(jRange)
    jRange = 1+annulus:N-1;
  else
    jRange = annulus+jRange(1:end-1);
  endif
  if isempty(kRange)
    kRange = 1:P-1;
  else
    kRange = kRange(1:end-1);
  endif
  M = numel(jRange)*numel(kRange);
  v0 = reshape(idx(i,  jRange,   kRange ),1,M);
  v1 = reshape(idx(i,  jRange,  kRange+1),1,M);
  v2 = reshape(idx(i, jRange+1, kRange+1),1,M);
  v3 = reshape(idx(i, jRange+1,  kRange ),1,M);
  if flip
    fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v3; v2; v1]);
  else
    fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v1; v2; v3]);
  endif
endfunction


%% Write IK boundary of a 3D grid (quad elements)
function writeIKboundary(fid,idx,j,iRange=[],kRange=[],flip=false)
  [N,~,P] = size(idx);
  if isempty(iRange)
    iRange = 1:N-1;
  else
    iRange = iRange(1:end-1);
  endif
  if isempty(kRange)
    kRange = 1:P-1;
  else
    kRange = kRange(1:end-1);
  endif
  M = numel(iRange)*numel(kRange);
  v0 = reshape(idx( iRange,  j,  kRange ),1,M);
  v1 = reshape(idx( iRange,  j, kRange+1),1,M);
  v2 = reshape(idx(iRange+1, j, kRange+1),1,M);
  v3 = reshape(idx(iRange+1, j,  kRange ),1,M);
  if flip
    fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v3; v2; v1]);
  else
    fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v1; v2; v3]);
  endif
endfunction


%% Write boundary of triangles (section of wedge block)
function writeTriPatch(fid,tri,idx,k,flip=false)
  v0 = idx(tri(:,1),k)';
  v1 = idx(tri(:,2),k)';
  v2 = idx(tri(:,3),k)';
  if flip
    fprintf(fid,"5  %i  %i  %i\n",[v0; v1; v2]);
  else
    fprintf(fid,"5  %i  %i  %i\n",[v0; v2; v1]);
  endif
endfunction


%% Write the top and bottom of the wedge block (quad elements)
function writeWedgeEnds(fid,idx,mask,kRange=[],bottom=true,top=true,flip=false)
  [N,P] = size(idx);
  if isempty(kRange)
    kRange = 1:P-1;
  else
    kRange = kRange(1:end-1);
  endif
  i0 = 1;  i3 = 0;
  while i3 != N
    for i2 = i0+1:N
      if mask(i2) != 0
        i1 = i2+1;
        for i3 = i1+1:N
          if mask(i3) != 0
            break;
          endif
        endfor
        break;
      endif
    endfor
    if bottom
      v0 = idx(i0, kRange );
      v1 = idx(i0,kRange+1);
      v2 = idx(i1,kRange+1);
      v3 = idx(i1, kRange );
      if flip
        fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v3; v2; v1]);
      else
        fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v1; v2; v3]);
      endif
    endif
    if top
      v0 = idx(i2, kRange );
      v1 = idx(i2,kRange+1);
      v2 = idx(i3,kRange+1);
      v3 = idx(i3, kRange );
      if flip
        fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v3; v2; v1]);
      else
        fprintf(fid,"9  %i  %i  %i  %i\n",[v0; v1; v2; v3]);
      endif
    endif
    i0 = i1;
  endwhile
endfunction
