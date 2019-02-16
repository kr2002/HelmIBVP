%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function setup boundary conditions for the 
% adjoint Helmholtz problem.
%
% The adjoint Helmholtz model:
%
% \Delta w + k^2(1+n(x)) w = 0  in \Omega
% w = g,  on \partial \Omega
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qmatrix,gmatrix,hmatrix,rmatrix] = HelmholtzBCAdj(p,e,u,time,g,ks,srcdetpair)

ne = size(e,2); % number of edges
qmatrix = zeros(1,ne);
gmatrix = zeros(1,ne);
hmatrix = ones(1,2*ne);
rmatrix = zeros(1,2*ne);

for k = 1:ne

	x1 = p(1,e(1,k)); % x at first point in segment
	y1 = p(2,e(1,k)); % y at first point in segment
	x2 = p(1,e(2,k)); % x at second point in segment
	y2 = p(2,e(2,k)); % y at second point in segment
	%xm = (x1 + x2)/2; % x at segment midpoint
	%ym = (y1 + y2)/2; % y at segment midpoint

    rmatrix(k)=g(k)*srcdetpair(ks,k);

end
for k=1:ne-1
    rmatrix(k+ne)=g(k+1)*srcdetpair(ks,k+1);
end
rmatrix(ne+ne)=g(1)*srcdetpair(ks,1);