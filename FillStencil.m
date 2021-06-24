%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [STC_A,bx,by]=FillStencil( ...
          a,c,NN,M,LIN,COL,index,X,Y,Xeta,Yeta,ALFA,BETA,GAMA,J)
%FILLSTENCIL populate matrix and source term
STC_A=zeros(3,3,NN);
bx=zeros(NN,1);   by=bx;
for n=1:NN
    % row and column
    i=LIN(n);   j=COL(n);
        
    % matrix (apply finite difference stencil)
    STC_A(:,:,n)=[-BETA(i,j)/2       GAMA(i,j)        BETA(i,j)/2; ...
                   ALFA(i,j) -2*(ALFA(i,j)+GAMA(i,j)) ALFA(i,j); ...
                   BETA(i,j)/2       GAMA(i,j)       -BETA(i,j)/2];
              
    % source term
    Q=a*exp(-c*(i-1));
    if i==M-1
        bx(n)=0.5*BETA(i,j)*(X(i+1,j+1)-X(i+1,j-1))-GAMA(i,j)*X(i+1,j)+ ...
              Q*J(i,j)^2*Xeta(i,j);
        by(n)=0.5*BETA(i,j)*(Y(i+1,j+1)-Y(i+1,j-1))-GAMA(i,j)*Y(i+1,j)+ ...
              Q*J(i,j)^2*Yeta(i,j);
    elseif i==2
        bx(n)=0.5*BETA(i,j)*(X(i-1,j-1)-X(i-1,j+1))-GAMA(i,j)*X(i-1,j)+ ...
              Q*J(i,j)^2*Xeta(i,j);
        by(n)=0.5*BETA(i,j)*(Y(i-1,j-1)-Y(i-1,j+1))-GAMA(i,j)*Y(i-1,j)+ ...
              Q*J(i,j)^2*Yeta(i,j);
    else
        bx(n)=Q*J(i,j)^2*Xeta(i,j);
        by(n)=Q*J(i,j)^2*Yeta(i,j);
    end
end
STC_A=reshape(STC_A,1,9*NN);
STC_A=STC_A(index);
end
