%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [X,Y,M] = Extrusion(Xin,Yin,data)

N = data.N;  R0 = data.L;  Rfar = data.targetFar;

Rprev = mean(sqrt(Xin(end-1,:).^2+Yin(end-1,:).^2))/data.chord;

% determine the radial locations
R = [Rprev R0];
while (Rfar-R(end))/(R(end)-R(end-1)) > 0.5
  i = numel(R);
  % constant rate
  dR = data.rate*(R(i)-R(i-1));
  % limit aspect ratio
  dT = 2*pi*R(i)/(N-1);
  dR = min(dR,dT*data.farAR);
  % append to R
  R = [R R(end)+dR];
endwhile
R(end) = Rfar;
R(end-1) = 0.5*(Rfar+R(end-2));
R = data.chord*R(3:end)';

% create circles
teta = atan2(Yin(end,:),Xin(end,:));

X = [Xin; R*cos(teta)];
Y = [Yin; R*sin(teta)];

% blend the transition
[M,~] = size(Xin);
for i=1:3
  for j=0:i
    for k=unique([M-j M+j])
      X(k,:) = 0.5*(X(k-1,:)+X(k+1,:));
      Y(k,:) = 0.5*(Y(k-1,:)+Y(k+1,:));
    endfor
  endfor
endfor
  
[M,~] = size(X);

endfunction
