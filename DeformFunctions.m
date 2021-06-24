%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md

1+1; clear ans % little hack to make this a "header" file

% Solve Ucp = M\Upt where M is the Wendland C2 RBF kernel built for xyz
function Ucp = SolveRBF(x,y,z,Upt)
  % dist = sqrt((x-x').^2+(y-y').^2+(z-z').^2);
  % M = max(1-dist,0).^4.*(1+4*dist);
  % Ucp = M\Upt;

  nPt = rows(x);
  nRhs = columns(Upt);

  fid = fopen("__coords.dat","w");
  fwrite(fid,[x y z],"double");
  fclose(fid);

  fid = fopen("__displs.dat","w");
  fwrite(fid,Upt,"double");
  fclose(fid);

  command = sprintf("./SolveRBF %d %d %d",nPt,nRhs,nproc);
  assert(!system(command),"External RBF solver failed.");

  fid = fopen("__ctrlpt.dat","r");
  Ucp = fread(fid,[nPt nRhs],"double");
  fclose(fid);

  system("rm __coords.dat  __displs.dat  __ctrlpt.dat");
endfunction

% Apply RBF deformation defined by Xcp and Ucp to points X
function [x,y,z] = ApplyRBF(Xcp,Ucp,X)
  nPt = rows(X);
  nCP = rows(Xcp);
  nDim = columns(Ucp);

  fid = fopen("__coords.dat","w");
  fwrite(fid,Xcp,"double");
  fclose(fid);

  fid = fopen("__ctrlpt.dat","w");
  fwrite(fid,Ucp,"double");
  fclose(fid);

  fid = fopen("__points.dat","w");
  fwrite(fid,X,"double");
  fclose(fid);

  command = sprintf("./ApplyRBF %d %d %d %d",nPt,nCP,nDim,nproc);
  assert(!system(command),"External RBF solver failed.");

  fid = fopen("__points.dat","r");
  x = fread(fid,nPt,"double");
  y = fread(fid,nPt,"double");
  z = fread(fid,nPt,"double");
  fclose(fid);

  system("rm __coords.dat __points.dat __ctrlpt.dat");
endfunction

% Apply RBF to all coordinates, xyz and uvw are column vectors
defXYZ = @(x,y,z,rbf)ApplyRBF([rbf.xcp  rbf.ycp  rbf.zcp]/rbf.radius,...
                              [rbf.ucp  rbf.vcp  rbf.wcp]/rbf.radius,...
                              [x y z]/rbf.radius);

% Apply RBF to only X and Y
defXY = @(x,y,z,rbf)ApplyRBF([rbf.xcp  rbf.ycp  rbf.zcp]/rbf.radius,...
                             [rbf.ucp  rbf.vcp]/rbf.radius,...
                             [x y z]/rbf.radius);

% Apply the deformation defined by fun and rbf to k XY slices of hex block XYZ
function [X,Y,Z] = defHexBlock(X,Y,Z,fun,rbf,k)
  [m,n,p] = size(X(:,:,1:k));
  x = reshape(X(:,:,1:k),m*n*p,1);
  y = reshape(Y(:,:,1:k),m*n*p,1);
  z = reshape(Z(:,:,1:k),m*n*p,1);
  [x,y,z] = fun(x,y,z,rbf);
  X(:,:,1:k) = reshape(x,m,n,p)*rbf.radius;
  Y(:,:,1:k) = reshape(y,m,n,p)*rbf.radius;
  Z(:,:,1:k) = reshape(z,m,n,p)*rbf.radius;
endfunction

% Same but for blocks of prisms, fun should be defXY or defXYZ
function [X,Y,Z] = defPrismBlock(X,Y,Z,fun,rbf,k)
  [n,p] = size(X(:,1:k));
  x = reshape(X(:,1:k),n*p,1);
  y = reshape(Y(:,1:k),n*p,1);
  z = reshape(Z(:,1:k),n*p,1);
  [x,y,z] = fun(x,y,z,rbf);
  X(:,1:k) = reshape(x,n,p)*rbf.radius;
  Y(:,1:k) = reshape(y,n,p)*rbf.radius;
  Z(:,1:k) = reshape(z,n,p)*rbf.radius;
endfunction
