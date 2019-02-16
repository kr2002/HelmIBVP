%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the value of the characteristic 
% function of the circle centered at circ(1:2) whose  
% radius is circ(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=ind_circ(x,circ)
N=length(x(1,:));
f(1:N)=0.0;
for k=1:N
    if norm(circ(1:2)'-x(:,k))<=circ(3)
        f(k)=1;
    end
end
f=f';  