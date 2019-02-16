%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helmholtz: Reconstruction algorithm for an inverse 
%            coefficient problem to the Helmholtz 
%            equation.
%            
%            The algorithm is based on least-square 
%            minimization.
%
% Author:    Kui Ren
% Address:   APAM, Columbia University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Helmholtz model:
%
% \Delta u + k^2(1+n(x)) u = 0  in \Omega
% u = f,  on \partial \Omega
%
% The measurement quantity:
% 
% g=\nu\cdot\nabla u  on \partial\Omega
%
% The data:
%
% (f_j, g_j), 1\le j\le N_s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The algorithm: minimize the functional
%
% \Phi:= 1/2*\sum_{j=1}^{N_s} \int_{\partial\Omega} (g_j-g_j^*)^2 ds
%
% Gradient of \Phi is computed with the adjoint state
% method where the adjoint problems are:
%
% \Delta w_j+ k^2(1+n)w_j=S in \Omega
% w_j = -r_j:=-(g_j-g_j^*)
%
% The gradient in direction \delta n is given as:
% \Phi' = \int_\Omega u_j w_j \delta n dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: 
%
% Step 1: Generate geometry using PDETOOL, and save the data
%
% Step 2: Setup simulations parameters accordingly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attention:
%
% We are using Dirichlet-to-Neumann data. Therefore 
% the code will evaluate gradient of solutions at the 
% boundary. The error of Neumann data at the boundary 
% could be very large if the simulations are not setup 
% careful enough. This could give very wrong simulation
% results.
%
% 1. Mesh should be fine enough to resolve the waves. 
%    When increasing k, we need to refine mesh as well.
%
% 2. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;  clear all; 

tic;
tb=toc;

% Load information on the domain
disp(' ');
disp(' ');
disp('Setting simulation geometry and parameters .......');
disp(' ');

load 'geo-1b1';

MaxIT=200; % max # of iterations allowed for the optimization algorithm
Ns=36; % number of illumination sources used

% Set wave number k: the size of domain is of order 1; we need to make 
%     sure that k is small enough so that there are enough grid points 
%     within each wavelength. Therefore, if we change k->2k, we need to
%     change mesh size h->h/4 at least.
k=4; % wave number (be very careful when increasing k, as noted above

% Create a Cartesian grid for inversion
dx=0.0125; x=0:dx:1;
dy=0.0125; y=0:dy:1;
Nx=length(x);
Ny=length(y);
% [X,Y]=meshgrid(x,y);

% Generate regular finite element mesh on rectangular geometry
% The mesh is not necessarily uniform. One can use other mesh here
[P,E,T]=poimesh(geo,Nx-1,Ny-1);
M=Nx*Ny; % total number of nodes in the mesh

srcinfo=SetSources(Ns);
[detinfo srcdetpair]=SetDetectors(P,E,Ns,srcinfo);
Nd=length(E(1,:)); % # of possible detectors, i.e number of boundary nodes, 
%                    also # of boundary edges

% Set true refractive index
reft=zeros(M,1);

%rec1=[0.2 0.4; 0.4 0.9];
%rec2=[1.4 1.7; 1.4 1.7];
%rec3=[1.0 1.8; 0.3 0.6];
circ1=[0.25 0.25 0.25];
circ2=[0.75 0.75 0.25];
%circ3=[0.5 1.2 0.3];
%circ4=[1.4 1.6 0.2];
%circ5=[1.7 0.4 0.2];
%circ6=[0.5 0.4 0.2];

%reft=0.4+0.2*ind_rec(P,rec1)+0.2*ind_rec(P,rec2)+0.2*ind_rec(P,rec3);
r1=(P(1,:)-0.25).^2+(P(2,:)-0.25).^2;
r1=sqrt(r1);
r2=(P(1,:)-0.75).^2+(P(2,:)-0.75).^2;
r2=sqrt(r2);
reft=0.5+0.5*cos(pi*r1/0.5)'.*ind_circ(P,circ1)+0.5*cos(pi*r2/0.5)'.*ind_circ(P,circ2);

% Interpolate to Cartesian grid
reftg=tri2grid(P,T,reft,x,y); % true value of n on Cartesian grid

figure;
pcolor(x,y,reftg); axis tight; colorbar('SouthOutside');
axis square; axis off; shading interp;
title('true refractive index n');
drawnow;

disp('Finished setting simulation geometry and parameters .......');

% Generating synthetic data
disp(' ');
disp(' ');
disp('Generating synthetic data .......');
disp(' ');

zerosrc=zeros(M,1); % interior source term is 0
measn=zeros(Nd,Ns); % measured noisy data
for ks=1:Ns
    
    % Solve the Helmholtz equation
    [ut meast]=HelmholtzSolve(P,E,T,k,reft,zerosrc,ks,Nd,srcinfo,detinfo);
       
    %utg=tri2grid(P,T,ut,x,y);
    %figure;
    %pcolor(x,y,utg); axis tight; colorbar('SouthOutside');
    %axis square; axis off; shading interp;
    %drawnow;
    %pause
    
    % Add multiplicative noise to data
    noiselevel=0.0;
    noise=noiselevel*2*(rand(1,Nd)-0.5);
	measn(:,ks)=meast.*(1+noise);
    
    %plot(meast);
    %pause;
    
    disp(['Synthetic data generated for source #: ' num2str(ks)]);
    disp('  ');

end
disp('Finished generating synthetic data .......');

% Setup initial guess
disp(' ');
disp(' ');
disp('Setting initial guess .......');
disp(' ');

ref0=0.5*ones(M,1); % initial guess of refractive index
ref0g=tri2grid(P,T,ref0,x,y); % interpolate onto Cartesian grid

figure;
pcolor(x,y,ref0g); axis tight; colorbar('SouthOutside');
axis square; axis off; shading interp;
%caxis([0.05 0.25]);
title('initial guess of refractive index n');
drawnow;

X0=ref0;

disp('Finished setting initial guess .......');

% This short part is only for debugging
%[f0 g0]=HelmholtzObj(X0,x,y,dx,dy,Nx,Ny,P,E,T,k,Ns,Nd,...
%                                BdaryNode,NormVecNode,ds,measn)                        
%g0g=tri2grid(P,T,g0,x,y);
%figure;
%pcolor(x,y,g0g); axis tight; colorbar('SouthOutside');
%axis square; shading interp;
%title('Gradient at initial guess');
%drawnow;

OptimMethod='UNCON'; % use uncontstrained (UNCON) or constrained (CON) minimization

% Setup the minimization algorithm
disp(' ');
disp(' ');
disp('Minimizing objective function .......');
disp(' ');

f=@(X) HelmholtzObj(X,x,y,dx,dy,Nx,Ny,P,E,T,k,Ns,Nd,srcinfo,detinfo,srcdetpair,measn);

if strcmp(OptimMethod,'UNCON')
    options=optimoptions(@fminunc,'Algorithm','quasi-newton', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    [X,fval,exitflag,output,grad]=fminunc(f,X0,options);
else
    % set inequality constraint
    Aieq=zeros(1,M);
    Bieq=0;
    % set equality constraint
    Aeq=zeros(1,M);
    Beq=0;
    % set upper and lower bounds
    LB=0.35*ones(1,M);
    UB=1.0*ones(1,M);

    options=optimoptions(@fmincon,'Algorithm','sqp', ...
    'Display','iter-detailed','GradObj','on','TolFun',1e-12,...
    'MaxIter',MaxIT);
    %options=optimset('Display','iter-detailed','GradObj','on','TolFun',1e-12,'MaxIter',MaxIT);
    %options = optimset('algorithm','sqp','maxfunevals',5000,'maxiter',100);
    %options = optimset(options,'tolx',1e-9,'tolcon',1e-9,'tolfun',1e-6);
    %options = optimset(options,'GradObj','on','GradConstr','off');
    
    [X,fval,exitflag,output,lambda]=fmincon(f,X0,Aieq,Bieq,Aeq,Beq,LB,UB,[],options);
    %[X,fval,exitflag,output]=fmincon(f,X0,zeros(M,M),zeros(M,1),[],[],LB,UB);
end
refr=X; % reconstructed refractive index

disp(' ');
disp(' ');
disp('Finished minimizing objective function .......');

disp(' ');
disp(' ');
disp('Plotting final results .......');
disp(' ');

% Plot reconstruction results
refrg=tri2grid(P,T,refr,x,y); % interpolate onto Cartesian grid
figure;
pcolor(x,y,refrg); axis tight; colorbar('SouthOutside');
axis square; axis off; shading interp;
title('reconstructed refractive index');
drawnow;

disp('Finished plotting final results .......');

% Save simulation results
save Exp01-info geo P E T srcinfo detinfo k MaxIT noiselevel dx dy
save Exp01-result reft ref0 refr -ASCII

te=toc;
disp(' ');
disp(' ');
disp(['The code run for: ' num2str(te-tb) ' seconds']);
disp(' ');
disp(' ');

% This last line is used to close MATLAB after the computation. It is 
% only used when runing the code in background.

%exit; % to exit MATLAB 