%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function chooses detectors to be used for each source
%
% Each boundary node is a detector. However, we choose 
% detectors that are not: (i) at the corners of the domain, 
% or (ii) too close to the source used.
% 
% The matrix detinfo is arranged as follows:
%
% detinfo(1,kd): the boundary node where detector is located
% detinfo(2,kd): the bdounary segment on which detector is located
% detinfo(3,kd): x conmponent of outer normal vector
% detinfo(4,kd): y conmponent of outer normal vector
% detinfo(5,kd): length of the boundary edge where detector is located
%
% The matrix srcdetpair is arranged as follows:
%
% srcdetpair(ks,kd) =1 if detector kd is used for source ks
% srcdetpair(ks,kd) =0 if detector kd is used for source ks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [detinfo srcdetpair]=SetDetectors(P,E,Ns,srcinfo)

Nd=length(E(1,:)); % total number of possible detectors (i.e boundary nodes)
detinfo=zeros(5,Nd);

detinfo(1,:)=E(1,:); % all boundary nodes can be detectors

% Find which segment of boundary each detector belongs to
for kd=1:Nd
    if abs(P(1,detinfo(1,kd))-0)<1e-8 % left
        if abs(P(2,detinfo(1,kd))-0)<1e-8
            detinfo(2,kd)=14;
        elseif abs(P(2,detinfo(1,kd))-1)<1e-8
            detinfo(2,kd)=3;
        else
            detinfo(2,kd)=4;
        end
    elseif abs(P(1,detinfo(1,kd))-1)<1e-8 % right
        if abs(P(2,detinfo(1,kd))-0)<1e-8
            detinfo(2,kd)=12;
        elseif abs(P(2,detinfo(1,kd))-1)<1e-8
            detinfo(2,kd)=3;
        else
            detinfo(2,kd)=2;
        end
    else
        if abs(P(2,detinfo(1,kd))-0)<1e-8 % bottom & top
            detinfo(2,kd)=1;
        else
            detinfo(2,kd)=3;
        end
    end
end  

% Find the outer normal vector and the length of each boundary edge 
% The boundary edges given in E are assumed to be ordered in the usual 
% convention: when traveling on the boundary, the domain is always on 
% the left. The outer normal vectors are found by rotating the tangent 
% vector by 90 degrees clockwise. This is done by applying the rotation 
% matrix:
% \cos\theta -\sin\theta
% \sin\theta \cos\theta
% with \theta=-90 since the rotation matrix rotates in the anticlock direction.
BdryEdgeVec=zeros(2,Nd);
for j=1:Nd
    ind1=E(1,j);
    ind2=E(2,j);
    vecx=P(1,ind2)-P(1,ind1);
    vecy=P(2,ind2)-P(2,ind1);
    veclength=sqrt(vecx^2+vecy^2);
    vecx=vecx/veclength;
    vecy=vecy/veclength;
    BdryEdgeVec(1,j)=vecy;
    BdryEdgeVec(2,j)=-vecx;
    detinfo(5,j)=veclength;
end

% Find the outer normal vectors of boundary nodes
% The normal vector at a node is defined here as the average of the outer 
% normal vector of the boundary edges joined together by the node
detinfo(3:4,1)=(BdryEdgeVec(:,Nd)+BdryEdgeVec(:,1))/2;
for j=2:Nd
    detinfo(3:4,j)=(BdryEdgeVec(:,j-1)+BdryEdgeVec(:,j))/2;
end

srcdetpair=ones(Ns,Nd);
% Exclude corner nodes
for kd=1:Nd
    if abs(P(1,detinfo(1,kd))-0)<1e-8 & abs(P(2,detinfo(1,kd))-0)<1e-8
        srcdetpair(:,kd)=0; % take out the (0,0) corner from detector list
    end
    if abs(P(1,detinfo(1,kd))-1)<1e-8 & abs(P(2,detinfo(1,kd))-0)<1e-8
        srcdetpair(:,kd)=0; % take out the (2,0) corner from detector list
    end
    if abs(P(1,detinfo(1,kd))-1)<1e-8 & abs(P(2,detinfo(1,kd))-1)<1e-8
        srcdetpair(:,kd)=0; % take out the (2,2) corner from detector list
    end
    if abs(P(1,detinfo(1,kd))-0)<1e-8 & abs(P(2,detinfo(1,kd))-1)<1e-8
        srcdetpair(:,kd)=0; % take out the (0,2) corner from detector list
    end
end
% Exclude nodes that are too close to the source location        
for ks=1:Ns
    for kd=1:Nd
        r=(srcinfo(1,ks)-P(1,detinfo(1,kd)))^2+(srcinfo(2,ks)-P(2,detinfo(1,kd)))^2;
        r=sqrt(r); % 
        if r<=0.2
            srcdetpair(ks,kd)=0; % take out the detecotrs too close to the source
        end
    end
end