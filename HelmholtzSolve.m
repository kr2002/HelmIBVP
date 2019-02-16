%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HelmholtzSolve: FEM solver for interior Helmholtz
%
% The Helmholtz model (assuming no resonance):
%
% \Delta u + k^2(1+n(x)) u = 0  in \Omega
% u = g,  on \partial \Omega
%
% The measurement quantity:
% 
% h=\nu\cdot\nabla u  on \partial\Omega
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u meas]=HelmholtzSolve(P,E,T,k,ref,f,...
    ks,Nd,srcinfo,detinfo)

% interpolation to triangle middle point
refm=pdeintrp(P,T,ref);
fm=pdeintrp(P,T,f);

% construct boundary conditions
pdebound =@(p,e,u,time)HelmholtzBC(p,e,[],[],ks,srcinfo,detinfo);
[Q,G,H,R] = assemb(pdebound,P,E);

% construct mass matrices
[K,M,F]=assema(P,T,-1,k^2*(1+refm),fm);

% solve the PDE
u = assempde(K,M,F,Q,G,H,R);
    
% solve PDEs (old way of solving PDE)
%u=assempde(geobc,P,E,T,-1,k^2*(1+refm),fm);

% compute measured data
[ux1, uy1]=pdegrad(P,T,u); % evaluate gradient of solution u
ux=pdeprtni(P,T,ux1); uy=pdeprtni(P,T,uy1); % intepolation to nodes
meas=zeros(1,Nd);
for j=1:Nd
    meas(j)=ux(detinfo(1,j))*detinfo(3,j)+uy(detinfo(1,j))*detinfo(4,j);
end