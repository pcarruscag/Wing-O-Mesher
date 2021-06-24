%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
% Visualization
if data.nDim == 2
  clf
  hold on
  for i=1:data.M
      plot(X(i,:),Y(i,:),'k')
  end

  for i=2:data.N
      plot(X(:,i),Y(:,i),'k')
  end
  axis equal
  hold off
  clear i
  return
endif

clf
hold on
mesh(squeeze(Xcyl(:,:,1)), squeeze(Ycyl(:,:,1)), squeeze(Zcyl(:,:,1)))
mesh(squeeze(Xcyl(1,:,:)), squeeze(Ycyl(1,:,:)), squeeze(Zcyl(1,:,:)))
mesh(Xq0(:,:,1), Yq0(:,:,1), Zq0(:,:,1))
mesh(Xq1(:,:,1), Yq1(:,:,1), Zq1(:,:,1))
trimesh(tri0, Xt0(:,1), Yt0(:,1), Zt0(:,1))
mesh(Xq2(:,:,1), Yq2(:,:,1), Zq2(:,:,1))
trimesh(tri1, Xt1(:,1), Yt1(:,1), Zt1(:,1))
mesh(Xq3(:,:,1), Yq3(:,:,1), Zq3(:,:,1))
trimesh(tri2, Xt2(:,1), Yt2(:,1), Zt2(:,1))
mesh(Xq4(:,:,1), Yq4(:,:,1), Zq4(:,:,1))
mesh(Xq5(:,:,1), Yq5(:,:,1), Zq5(:,:,1))
caxis([0 1e6]); colormap(gray);
axis equal

figure
hold on
mesh(squeeze(Xcyl(1,:,end-1:end)),...
     squeeze(Ycyl(1,:,end-1:end)),...
     squeeze(Zcyl(1,:,end-1:end)))
mesh(Xq0(:,:,1), Yq0(:,:,1), Zq0(:,:,1))
mesh(Xq1(:,:,1), Yq1(:,:,1), Zq1(:,:,1))
trimesh(tri0, Xt0(:,1), Yt0(:,1), Zt0(:,1))
mesh(Xq2(:,:,1), Yq2(:,:,1), Zq2(:,:,1))
trimesh(tri1, Xt1(:,1), Yt1(:,1), Zt1(:,1))
mesh(Xq3(:,:,1), Yq3(:,:,1), Zq3(:,:,1))
trimesh(tri2, Xt2(:,1), Yt2(:,1), Zt2(:,1))
mesh(Xq4(:,:,1), Yq4(:,:,1), Zq4(:,:,1))
mesh(Xq5(:,:,1), Yq5(:,:,1), Zq5(:,:,1))
caxis([0 1e6]); colormap(gray);
axis equal

figure
hold on
mesh(squeeze(Xan(:,:,1)),   squeeze(Yan(:,:,1)),   squeeze(Zan(:,:,1)))
mesh(squeeze(Xan(end,:,:)), squeeze(Yan(end,:,:)), squeeze(Zan(end,:,:)))
mesh(Xq0(:,:,end), Yq0(:,:,end), Zq0(:,:,end))
mesh(Xq1(:,:,end), Yq1(:,:,end), Zq1(:,:,end))
trimesh(tri0, Xt0(:,end), Yt0(:,end), Zt0(:,end))
mesh(Xq2(:,:,end), Yq2(:,:,end), Zq2(:,:,end))
trimesh(tri1, Xt1(:,end), Yt1(:,end), Zt1(:,end))
mesh(Xq3(:,:,end), Yq3(:,:,end), Zq3(:,:,end))
trimesh(tri2, Xt2(:,end), Yt2(:,end), Zt2(:,end))
mesh(Xq4(:,:,end), Yq4(:,:,end), Zq4(:,:,end))
mesh(Xq5(:,:,end), Yq5(:,:,end), Zq5(:,:,end))
caxis([0 1e6]); colormap(gray);
axis equal

figure
hold on
mesh(squeeze(Xan(:,2,:)),squeeze(Zan(:,2,:)),squeeze(Yan(:,2,:)))
mesh(squeeze(Xcyl(:,2,:)),squeeze(Zcyl(:,2,:)),squeeze(Ycyl(:,2,:)))
i = ceil(rows(Xq0)/2);
mesh(squeeze(Xq0(i,:,:)),squeeze(Zq0(i,:,:)),squeeze(Yq0(i,:,:)))
i = ceil(rows(Xq1)/2);
mesh(squeeze(Xq1(i,:,:)),squeeze(Zq1(i,:,:)),squeeze(Yq1(i,:,:)))
caxis([0 1e6]); colormap(gray);
axis equal
clear i
