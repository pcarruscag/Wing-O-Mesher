%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [x,y,z,Nwing] = GetWingSurfPoints(endOffset, Xcyl,Ycyl,Zcyl,...
  Xq0,Yq0,Zq0, Xq1,Yq1,Zq1, Xt0,Yt0,Zt0, Xq2,Yq2,Zq2, Xt1,Yt1,Zt1,...
  Xq3,Yq3,Zq3, Xt2,Yt2,Zt2, Xq4,Yq4,Zq4, Xq5,Yq5,Zq5)

N = numel(Xcyl(1,2:end-1,max(1,end-endOffset):end-1));
x = reshape(squeeze(Xcyl(1,2:end-1,max(1,end-endOffset):end-1)),N,1);
y = reshape(squeeze(Ycyl(1,2:end-1,max(1,end-endOffset):end-1)),N,1);
z = reshape(squeeze(Zcyl(1,2:end-1,max(1,end-endOffset):end-1)),N,1);
Nwing = N;

N = numel(Xq0(:,:,1));
x = [x; reshape(Xq0(:,:,1),N,1)];
y = [y; reshape(Yq0(:,:,1),N,1)];
z = [z; reshape(Zq0(:,:,1),N,1)];

N = numel(Xq1(:,2:end-1,1));
x = [x; reshape(Xq1(:,2:end-1,1),N,1)];
y = [y; reshape(Yq1(:,2:end-1,1),N,1)];
z = [z; reshape(Zq1(:,2:end-1,1),N,1)];

x = [x; Xt0(:,1)];
y = [y; Yt0(:,1)];
z = [z; Zt0(:,1)];

N = numel(Xq2(:,2:end-1,1));
x = [x; reshape(Xq2(:,2:end-1,1),N,1)];
y = [y; reshape(Yq2(:,2:end-1,1),N,1)];
z = [z; reshape(Zq2(:,2:end-1,1),N,1)];

x = [x; Xt1(:,1)];
y = [y; Yt1(:,1)];
z = [z; Zt1(:,1)];

N = numel(Xq3(:,2:end-1,1));
x = [x; reshape(Xq3(:,2:end-1,1),N,1)];
y = [y; reshape(Yq3(:,2:end-1,1),N,1)];
z = [z; reshape(Zq3(:,2:end-1,1),N,1)];

x = [x; Xt2(:,1)];
y = [y; Yt2(:,1)];
z = [z; Zt2(:,1)];

N = numel(Xq4(:,2:end-1,1));
x = [x; reshape(Xq4(:,2:end-1,1),N,1)];
y = [y; reshape(Yq4(:,2:end-1,1),N,1)];
z = [z; reshape(Zq4(:,2:end-1,1),N,1)];

N = numel(Xq5(:,:,1));
x = [x; reshape(Xq5(:,:,1),N,1)];
y = [y; reshape(Yq5(:,:,1),N,1)];
z = [z; reshape(Zq5(:,:,1),N,1)];

endfunction
