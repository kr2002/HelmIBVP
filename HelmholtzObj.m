%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The objective functional for the minimization problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f g]=HelmholtzObj(X,x,y,dx,dy,Nx,Ny,P,E,T,k,...
                            Ns,Nd,srcinfo,detinfo,srcdetpair,meas)

M=Nx*Ny;
refc=X; % current value of the refractive index

f=0.0;
g=zeros(M,1);

zerosrc=zeros(M,1);
for ks=1:Ns
    
    pred=zeros(1,Nd); % predicted data on measurement locations
    rz=zeros(1,Nd); % residual on measurement locations
 
    [uc pred]=HelmholtzSolve(P,E,T,k,refc,zerosrc,ks,Nd,srcinfo,detinfo);
        
    %ucg=tri2grid(P,T,uc,x,y);
    %figure;
    %pcolor(x,y,ucg); axis tight; colorbar('SouthOutside');
    %axis square; axis off; shading interp;
    %drawnow;
    
    measR=meas(:,ks);
    %rz=(pred-measR')./measR';
    rz=(pred-measR');
    % calculate the objective function (the part for source ks)
    f=f+0.5*sum(rz.^2.*detinfo(5,:).*srcdetpair(ks,:)); % rectangular rule for boundary integral
    
    if nargout > 1 % calculate gradient
        
        % solve the adjoint Helmholtz problem
        %src=-rz./measR';
        adjsrc=-rz;
        wc=HelmholtzSolveAdj(P,E,T,k,refc,adjsrc,ks,srcdetpair);
        % calculate the gradient w.r.t the refractive index
        g=g+k^2*uc.*wc*dx*dy;
 
    end
    
end

beta=1e-16; % the regularization parameter, change accoring to noise level

% Add regularization to the objective function
[Gx,Gy] = pdegrad(P,T,refc);
Gx1=pdeprtni(P,T,Gx); Gy1=pdeprtni(P,T,Gy);
f=f+0.5*beta*sum(Gx1.^2+Gy1.^2)*dx*dy;

% Add regularization term to the gradient
if nargout >1   
    % The Laplacian part
    [Gxx, Gxy]=pdegrad(P,T,Gx1); [Gyx, Gyy]=pdegrad(P,T,Gy1);
    Gx2=pdeprtni(P,T,Gxx); Gy2=pdeprtni(P,T,Gyy);
    DeltaGamma=Gx2+Gy2;
    g=g-beta*DeltaGamma*dx*dy;
    
    % The boundary integral part: this part should be removed if boundary 
    %      value of n is not to be reconstructed
    dndn=zeros(1,Nd);
    for j=1:Nd
        dndn(j)=Gx1(detinfo(1,j))*detinfo(3,j)+Gy1(detinfo(1,j))*detinfo(4,j);
    end
    g=g+beta*sum(dndn.*detinfo(5,:));

end