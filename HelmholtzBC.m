%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function setup boundary conditions for the 
% Helmholtz problem.
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
function [qmatrix,gmatrix,hmatrix,rmatrix] = HelmholtzBC(p,e,u,time,ks,srcinfo,detinfo)

ne = size(e,2); % number of edges on the domain boundary
qmatrix = zeros(1,ne);
gmatrix = qmatrix;
hmatrix = ones(1,2*ne);
rmatrix = zeros(1,2*ne);

% The sources we use are Gaussian sources (to mimic point sources)
xs=srcinfo(1,ks);
ys=srcinfo(2,ks);
srcseg=srcinfo(3,ks);
for k = 1:ne

	x1 = p(1,e(1,k)); % x at first point in segment
	y1 = p(2,e(1,k)); % y at first point in segment
	x2 = p(1,e(2,k)); % x at second point in segment
	y2 = p(2,e(2,k)); % y at second point in segment
	%xm = (x1 + x2)/2; % x at segment midpoint
	%ym = (y1 + y2)/2; % y at segment midpoint

    if detinfo(2,k)==srcseg % if the edge lives on the same side with the source
        rmatrix(k) = exp(-((x1-xs)^2+(y1-ys)^2)/0.01);
        rmatrix(k+ne) = exp(-((x2-xs)^2+(y2-ys)^2)/0.01);
    end

end
rmatrix=1.0e0*rmatrix; % change strength of sources