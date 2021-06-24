%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [X,Y,N] = Wake(Xin,Yin,data)

Nin = data.N;  M = data.M;

if data.wake > 1
  % 5 new grid points will be created
  N = Nin+5;
  X = zeros(M,N+1); Y = zeros(M,N+1);

  % copy the original points
  X(:,1:Nin/2) = Xin(:,1:Nin/2);  Y(:,1:Nin/2) = Yin(:,1:Nin/2);

  X(:,Nin/2+2) = Xin(:,Nin/2+1);  Y(:,Nin/2+2) = Yin(:,Nin/2+1);
  X(:,Nin/2+6) = Xin(:,Nin/2+2);  Y(:,Nin/2+6) = Yin(:,Nin/2+2);

  X(:,(Nin/2+8):end) = Xin(:,(Nin/2+3):end);
  Y(:,(Nin/2+8):end) = Yin(:,(Nin/2+3):end);

  % interpolate new points
  % split near trailing edge
  X(:,Nin/2+1) = 0.62*X(:,Nin/2+2)+0.38*X(:,Nin/2);
  Y(:,Nin/2+1) = 0.62*Y(:,Nin/2+2)+0.38*Y(:,Nin/2);

  X(:,Nin/2+7) = 0.62*X(:,Nin/2+6)+0.38*X(:,Nin/2+8);
  Y(:,Nin/2+7) = 0.62*Y(:,Nin/2+6)+0.38*Y(:,Nin/2+8);

  % split trailing edge
  X(:,Nin/2+4) = 0.5*(X(:,Nin/2+2)+X(:,Nin/2+6));
  Y(:,Nin/2+4) = 0.5*(Y(:,Nin/2+2)+Y(:,Nin/2+6));

  X(:,Nin/2+3) = 0.62*X(:,Nin/2+4)+0.38*X(:,Nin/2+2);
  Y(:,Nin/2+3) = 0.62*Y(:,Nin/2+4)+0.38*Y(:,Nin/2+2);

  X(:,Nin/2+5) = 0.62*X(:,Nin/2+4)+0.38*X(:,Nin/2+6);
  Y(:,Nin/2+5) = 0.62*Y(:,Nin/2+4)+0.38*Y(:,Nin/2+6);
else
  % only split TE
  N = Nin+1;
  X = zeros(M,N+1); Y = zeros(M,N+1);

  % copy the original points
  X(:,1:(Nin/2+1)) = Xin(:,(1:Nin/2+1));
  Y(:,1:(Nin/2+1)) = Yin(:,(1:Nin/2+1));

  X(:,(Nin/2+3):end) = Xin(:,(Nin/2+2):end);
  Y(:,(Nin/2+3):end) = Yin(:,(Nin/2+2):end);

  % interpolate new points
  % split trailing edge
  X(:,Nin/2+2) = 0.5*(X(:,Nin/2+1)+X(:,Nin/2+3));
  Y(:,Nin/2+2) = 0.5*(Y(:,Nin/2+1)+Y(:,Nin/2+3));
endif

endfunction
