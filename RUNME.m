%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
close all
clear
clc

%% GLOBAL PARAMETERS %%

INPUT_FILE = "mesh_settings";
INPUT_TYPE = "TXT";  % get input from text (TXT, .m) or load binary (BIN, .dat)
EXPORT = false;  % write mesh files
RAM_LIMIT = 14;  % abort computation if more than XX GB would be used

%% LOAD SETTINGS %%

if strcmp(INPUT_TYPE,"TXT")
  str2func(INPUT_FILE)();
  %save("-binary",[INPUT_FILE ".dat"],"data");
elseif strcmp(INPUT_TYPE,"BIN")
  load([INPUT_FILE ".dat"]);
else
  assert(false,"Unknown input type.");
endif

%% MESH GENERATION %%

printf('### 2D MESH GENERATION ###\n'); fflush(stdout);
[X,Y,iters,residual] = MeshGenerator(data);
[X,Y,data.M,data.N] = Refinement(X,Y,data);
[X,Y,data.M,growthRate] = Inflation(X,Y,data);
if data.nDim == 2 && data.targetFar > data.L
  [X,Y,data.M] = Extrusion(X,Y,data);
endif
if data.wake > 0
  [X,Y,data.N] = Wake(X,Y,data);
endif
if data.nDim == 3
  printf('\n### 3D MESH GENERATION ###\n');
  printf('  Extruding 2D mesh............ '); fflush(stdout);
  tic;
  Cylinder;
  MeshTip;
  SplitSurfaces;
  printf('%.2fs\n',toc);
  printf('  Spanwise interface patches: %i\n', surfSplits.number);
  fflush(stdout);

  if data.smoothTip || data.sweep!=0 || data.taper!=1 || data.extDef
    printf('  3D RBF mesh deformation:\n'); fflush(stdout);
    DeformFunctions;
    SmoothCorner;
    Sweep;
  endif
endif

%% ANALYSIS %%

j = round(data.N/4);
aspectRatio2 = sqrt((X(1,j+1)-X(1,j))^2+(Y(1,j+1)-Y(1,j))^2)/data.firstLayer;
if data.nDim == 3
  aspectRatio3 = max(diff(Zcyl(1,1,:)))/data.firstLayer;
endif
clear j
PlotMesh;

%% EXPORT %%

if EXPORT
  printf('\n### FLUID MESH EXPORT ###\n'); fflush(stdout);
  tic;
  if data.nDim == 2
    Save2D;
  else
    SaveFunctions;
    Save3D;
  endif
  [~,size] = system('ls -s --block-size=1 mesh.su2');
  size = str2num(strtok(size,' '))/2^20;
  printf('  Wrote %.0fMB in %.2fs\n',size,toc); fflush(stdout);
  clear size
endif

%% SOLID MESH %%

if data.fsi
  printf('\n### %iD SOLID MESH ###\n',data.nDim); fflush(stdout);
  tic;
  if data.nDim == 2
    SolidDomain2D;
  else
    SolidDomain3D;
  endif
  if EXPORT
    assert(!system(['./MergeMeshes.py ', num2str(data.nDim)]),
           'Could not merge fluid and solid mesh.');
  endif
  printf('  Time: %.2fs\n',toc);
endif

printf('\nMesh generation SUCCESSFUL.\n\n');
