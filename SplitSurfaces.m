%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
%% Work out split locations along span to limit the number of interface points

% no need to split if not FSI
if !data.fsi
  surfSplits.number = 1;
  surfSplits.sections_f = [1 Nz];
  surfSplits.coords = squeeze(Zcyl(1,1,surfSplits.sections_f));
  return
endif

% in general the structural side is the limiting factor, generate solid z points

if data.match
  % matching interface, use z locations of fluid sections
  z = Zcyl(1,2,:);  z(end) = data.chord*data.Lz;
else
  % non-matching, determine z locations from target element aspect ratios
  szAve = data.chord*data.t*0.01/data.nTh1;
  szRoot = data.arRoot*szAve;
  szTip = data.arTip*szAve*data.taper;
  Lz = data.chord*data.Lz;
  t0 = 0.5*(szRoot+szTip)/Lz;
  a = (szRoot/Lz-t0^2)/(t0-t0^2);
  t = 0:1/round(1/t0):1;
  z = zeros(1,1,numel(t));
  z(1,1,:) = (a*t+(1-a)*t.^2)*Lz;

  clear szAve szRoot szTip t0 a t
endif

Nz_s = numel(z);

nSurfPts = (cut9-cut7+1)*Nz_s;
nPatch = ceil(nSurfPts/data.maxInterfacePts);

icuts = 1:ceil(Nz_s/nPatch):Nz_s;
if(icuts(end) != Nz_s) icuts = [icuts Nz_s]; endif
zcuts = z(icuts);

% cut last fluid patch in half if it goes over limit
dist = abs(squeeze(Zcyl(1,1,:)-zcuts(end-1)));
[~,ifluid] = min(dist);
if (cut9-cut7+1)*(Nz-ifluid+1) > data.maxInterfacePts
  znew = Zcyl(1,1,round((ifluid+Nz)/2));
  % find closest solid point to insert cut
  dist = abs(squeeze(z)-znew);
  [~,isolid] = min(dist);
  icuts = [icuts(1:end-1) isolid Nz_s];
  zcuts = z(icuts);
  nPatch += 1;
endif

% find the closest z point on the fluid size for each cut, use that as location
% instead and adjust the solid z distribution accordingly

% normalized distribution for each section and closest fluid point
icuts_f = [];

for i=1:nPatch
  i0 = icuts(i); i1 = icuts(i+1);
  z0 = zcuts(i); z1 = zcuts(i+1);
  z(i0:i1) = [0 (z(i0+1:i1)-z0)/(z1-z0)];

  dist = abs(squeeze(Zcyl(1,1,:)-z0));
  [~,ifluid] = min(dist);

  icuts_f = [icuts_f ifluid];
  zcuts(i) = Zcyl(1,1,ifluid);
endfor

icuts_f = [icuts_f Nz];

% re-dim distribution
for i=1:nPatch
  i0 = icuts(i); i1 = icuts(i+1);
  z0 = zcuts(i); z1 = zcuts(i+1);
  z(i0:i1) = [0 z(i0+1:i1-1) 1]*(z1-z0)+z0;
endfor

surfSplits.number = nPatch;
surfSplits.sections_f = icuts_f;
surfSplits.sections_s = icuts;
surfSplits.coords = zcuts;
surfSplits.zdist = z;

clear dist i0 i1 icuts icuts_f ifluid isolid nPatch nSurfPts z z0 z1 zcuts znew
