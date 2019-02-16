%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HelmholtzSolve: FEM solver for interior Helmholtz
%
% The Helmholtz model (assuming no resonance):
%
% \Delta u + k^2(1+n(x)) u = f  in \Omega
% u = g,  on \partial \Omega
%
% The measurement quantity:
% 
% h=\nu\cdot\nabla u  on \partial\Omega
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=HelmholtzSolveAdj(P,E,T,k,ref,g,ks,srcdetpair)

% interpolation to triangle middle point
refm=pdeintrp(P,T,ref);

% construct boundary conditions
pdebound =@(p,e,u,time)HelmholtzBCAdj(p,e,[],[],g,ks,srcdetpair);
[Q,G,H,R] = assemb(pdebound,P,E);

% construct mass matrices
[K,M,F]=assema(P,T,-1,k^2*(1+refm),0);

% solve the PDE
u = assempde(K,M,F,Q,G,H,R);
    
% solve PDEs (old way of solving the PDE)
%u=assempde(geobc,P,E,T,-1,k^2*(1+refm),fm);