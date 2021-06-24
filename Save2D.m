%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
%% SAVE %%

fid = fopen("mesh.su2","w");
fprintf(fid,"NDIME= 2\n");

% number of vertices and elements
Nv = (data.N-1)*data.M;
Ne = (data.N-1)*(data.M-1);

printf("  Writing fluid mesh...\n");
printf("   Num. vertices: %i\n",Nv);
printf("   Num. elements: %i\n",Ne);

fprintf(fid,"NELEM= %i\n",Ne);

for i=1:(data.M-1)
  for j=2:data.N
    % the element index
    eli = (i-1)*(data.N-1)+j-2;
    % the vertices that define it
    v0 = (i-1)*(data.N-1)+j-2;
    v3 = i*(data.N-1)+j-2;
    if j==data.N % wrap around
      v1 = (i-1)*(data.N-1);
      v2 = i*(data.N-1);
    else
      v1 = (i-1)*(data.N-1)+j-1;
      v2 = i*(data.N-1)+j-1;
    endif
    
    fprintf(fid,"9  %i  %i  %i  %i  %i\n",v0,v1,v2,v3,eli);
  endfor
endfor

fprintf(fid,"NPOIN= %i\n",Nv);

for i=1:data.M
  for j=2:data.N
    % the vertex index
    vi = (i-1)*(data.N-1)+j-2;
    fprintf(fid,"%e  %e  %i\n",X(i,j),Y(i,j),vi);
  endfor
endfor

if data.splitSurfaceAt < 1
  fprintf(fid,"%s","NMARK= 2\n");
else
  fprintf(fid,"%s","NMARK= 4\n");
endif

fprintf(fid,"MARKER_TAG= farfield\n");
fprintf(fid,"MARKER_ELEMS= %i\n",data.N-1);

for j=2:data.N
  v0 = (data.M-1)*(data.N-1)+j-2;
  if j==data.N
    v1 = (data.M-1)*(data.N-1);
  else
    v1 = (data.M-1)*(data.N-1)+j-1;
  endif
  
  fprintf(fid,"3  %i  %i\n",v0,v1);
endfor

if data.splitSurfaceAt < 1
  fprintf(fid,"MARKER_TAG= airfoil\n");
  fprintf(fid,"MARKER_ELEMS= %i\n",data.N-1);

  for j=2:data.N-1
    fprintf(fid,"3  %i  %i\n",j-2,j-1);
  endfor
  fprintf(fid,"3  %i  %i\n",data.N-2,0);
  
else
  fprintf(fid,"MARKER_TAG= leading_edge\n");
  fprintf(fid,"MARKER_ELEMS= %i\n",2*data.splitSurfaceAt-2);
  
  for j=2:data.splitSurfaceAt
    fprintf(fid,"3  %i  %i\n",j-2,j-1);
  endfor
  fprintf(fid,"3  %i  %i\n",data.N-2,0);
  for j=3:data.splitSurfaceAt
    fprintf(fid,"3  %i  %i\n",data.N-j,data.N-j+1);
  endfor
  
  tmp = (data.N-1)/2-data.splitSurfaceAt+1;
  fprintf(fid,"MARKER_TAG= suction_side\n");
  fprintf(fid,"MARKER_ELEMS= %i\n",tmp);
  
  for j=1:tmp
    fprintf(fid,"3  %i  %i\n",data.splitSurfaceAt+j-2,data.splitSurfaceAt+j-1);
  endfor
  
  fprintf(fid,"MARKER_TAG= pressure_side\n");
  fprintf(fid,"MARKER_ELEMS= %i\n",tmp);
  
  for j=1:tmp
    fprintf(fid,"3  %i  %i\n",data.splitSurfaceAt+tmp+j-2,...
                              data.splitSurfaceAt+tmp+j-1);
  endfor
endif

fclose(fid);

clear Ne Nv ans eli fid i j v0 v1 v2 v3 vi tmp
