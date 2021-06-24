%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [X,Y,M,N] = Refinement(Xin,Yin,data)
% initialize
margin = 2;
X = [Xin(:,end-margin:end-1)  Xin(:,2:end)  Xin(:,3:2+margin)];
Y = [Yin(:,end-margin:end-1)  Yin(:,2:end)  Yin(:,3:2+margin)];
[M,N] = size(X);

for k=1:data.refineFactor
  % new sizes
  M = 2*M-1;  N = 2*N;
  margin *= 2;
  
  % interpolate
  isrc = 1:2:M;  jsrc = [1:2:(N/2-1)  (N/2+2):2:N];
  X = interp2(jsrc',isrc,X,(1:N)',1:M,'spline');
  Y = interp2(jsrc',isrc,Y,(1:N)',1:M,'spline');
endfor

% remove margin
X = X(:,1+margin:end-margin);
Y = Y(:,1+margin:end-margin);
N -= 2*margin;

% recover periodicity
X = [X(:,N-1) X];  Y = [Y(:,N-1) Y];

endfunction
