%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function tri = Triangulate(x,y)

tri = delaunay(x,y);

% filter out skewed triangles (good ones have 1 side with consecutive nodes)
keep = sort(tri,2);
keep = diff(keep,1,2);
keep = or(keep(:,1) == 1, keep(:,2) == 1);
tri = tri(keep,:);

% check orientation
for i=1:rows(tri)
  o = [x(tri(i,1)) y(tri(i,1)) 0];
  a = [x(tri(i,2)) y(tri(i,2)) 0] - o;
  b = [x(tri(i,3)) y(tri(i,3)) 0] - o;
  if cross(a,b)(3) < 0
    tmp = tri(i,2);
    tri(i,2) = tri(i,3);
    tri(i,3) = tmp;
  endif
endfor

endfunction
