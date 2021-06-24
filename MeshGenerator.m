%  Copyright (C) 2018-2021  Pedro Gomes
%  See full notice in NOTICE.md
%
function [X,Y,iter,residual] = MeshGenerator(data)

N=data.N; M=data.M; a=data.a; c=data.c;
ri=data.ri; rf=data.rf; kf=data.kf;
tol=data.linSolTolerance; imax=data.linSolMaxIters;

printf("  Generating elliptic mesh..... "); fflush(stdout);
tic
%% Initialization
% Number of equations
NN=(N-1)*(M-2);

% Matrices and vectors for numbering
NUM=reshape(1:NN,N-1,M-2)';
NUM=[NUM(:,N-1) NUM NUM(:,1)];
NUM=[zeros(1,N+1); NUM; zeros(1,N+1)];
COL=repmat(1:(N-1),1,M-2)+1;
LIN=reshape(repmat((1:(M-2))',1,N-1)',1,NN)+1;    

% Vector of nodes and equations
STC_N=zeros(3,3,NN);    
for n=1:NN
    i=LIN(n);   j=COL(n);
    STC_N(:,:,n)=NUM((i-1):(i+1),(j-1):(j+1));
end
STC_N=reshape(STC_N,1,9*NN);
index=(STC_N~=0);
STC_N=STC_N(index);
EQ=reshape(repmat(1:NN,9,1),1,9*NN);
EQ=EQ(index);
    
% Initial algebraic mesh
% Polar coordinates of the circles
R=(0.1:1/(M-1):1.1)/1.1;
R=0.5+(data.L-0.5)*R.^2;
R=repmat(R',1,N);
teta=2*pi*(0:1/(N-1):1);
blend = (0:(M-1))'/(M-1)*data.wakeCluster;
teta += (1-1.5*data.wakeUpDown)*blend*sin((1+data.wakeUpDown)*teta);
COS=-cos(teta);
SIN=sin(teta);

% Initial coordinate matrices
X = COS.*R;
Yscl=(1+7*exp(-0.4*(1:M))).^-1; % the circles are squached near the airfoil
Y=SIN.*R.*repmat(Yscl',1,N);
if data.wakeUpDown
  k = 0;
else
  k = data.wakeCluster;
endif
[X(1,:),Y(1,:)]=NACA4(data.m,data.p,data.t,N,k,data.teCoarsening);
X=[X(:,N-1) X];     Y=[Y(:,N-1) Y];
x0=reshape(X(2:(M-1),2:N)',NN,1);
y0=reshape(Y(2:(M-1),2:N)',NN,1);

%% Iterations
for k=1:data.maxIters
    % Mesh metrics by central differences
    Xqsi=0.5*[zeros(M,1)    X(:,3:(N+1))-X(:,1:(N-1)) zeros(M,1)];
    Xeta=0.5*[zeros(1,N+1); X(3:M,:)-X(1:(M-2),:);    zeros(1,N+1)];
    Yqsi=0.5*[zeros(M,1)    Y(:,3:(N+1))-Y(:,1:(N-1)) zeros(M,1)];
    Yeta=0.5*[zeros(1,N+1); Y(3:M,:)-Y(1:(M-2),:);    zeros(1,N+1)];
    % First order near the trailing edge
    te1 = N/2+1;  te2 = te1+1;
    Xqsi(:,te1) = X(:,te1)-X(:,te1-1);
    Xqsi(:,te2) = X(:,te2+1)-X(:,te2);
    Yqsi(:,te1) = Y(:,te1)-Y(:,te1-1);
    Yqsi(:,te2) = Y(:,te2+1)-Y(:,te2);
    
    J=Xqsi.*Yeta-Xeta.*Yqsi;
    
    % Constants for the equations
    ALFA=Xeta.^2+Yeta.^2;
    BETA=Xqsi.*Xeta+Yqsi.*Yeta;
    GAMA=Xqsi.^2+Yqsi.^2;
    
    % Relaxation factor
    rel=ri+(rf-ri)/kf*k;
    if rel>rf
        rel=rf;
    end
    
    % Populate matrix and source term
    [STC_A,bx,by]=FillStencil( ...
                  a,c,NN,M,LIN,COL,index,X,Y,Xeta,Yeta,ALFA,BETA,GAMA,J);
    A=sparse(EQ,STC_N,STC_A,NN,NN,9*NN);
    
    % New coordinates
    if k==1
        x1=A\bx;
        y1=A\by;
    else
        [iluL,iluU] = ilu(A);
        [x1,~]=bicgstab(A,bx,tol,imax,iluL,iluU,x0);
        [y1,~]=bicgstab(A,by,tol,imax,iluL,iluU,y0);
    end

    Xnext=reshape(x1,N-1,M-2)';
    Xnext=cat(2,Xnext(:,N-1),Xnext,Xnext(:,1));
    Xnext=cat(1,X(1,:),Xnext,X(M,:));
    X=rel*Xnext+(1-rel)*X;

    Ynext=reshape(y1,N-1,M-2)';
    Ynext=cat(2,Ynext(:,N-1),Ynext,Ynext(:,1));
    Ynext=cat(1,Y(1,:),Ynext,Y(M,:));
    Y=rel*Ynext+(1-rel)*Y;
    
    % Convergence criteria
    h=max(abs(x1-x0))+max(abs(y1-y0));
    if h<data.tolerance;
        break
    end
    y0=y1;  x0=x1;
end
X=X*data.chord;  Y=Y*data.chord;

iter=k;  residual=h;

printf('%.2fs\n',toc); fflush(stdout);

endfunction
