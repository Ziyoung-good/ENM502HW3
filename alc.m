function [u,lambda,deltas]=alc(u,lambda,u0,lambda0)

deltas=0.5;                               
dx=1/29;                     
dy=1/29;                     
x=0:dx:1;                        
y=0:dy:1; 
J=zeros(30*30+1);
b=zeros(30*30+1,1);
dRdl=zeros(30*30,1);
dEdu=zeros(30*30,1);
    
%compute J 
for i=1:30
    for j=1:30
        index=(j-1)*(30)+i;
        dEdu(index)=-2*(u(index)-u0(index));
        if ((0<y(j))&&(y(j)<1)&&(0<x(i))&&(x(i)<1))
            dRdl(index)=u(index)*(1+u(index));
            J(index,index)=-2/(dx*dx)-2/(dy*dy)+lambda*(1+2*u(index));
            J(index,index-1)=1/dx/dx;
            J(index,index+1)=1/dx/dx;
            J(index,index+30)=1/dy/dy;
            J(index,index-30)=1/dy/dy;
        else
            dRdl(index)=0;
            J(index,index)=1;
        end
    end
end
% compute the J
J(30*30+1,1:30*30)=dEdu;
J(1:30*30,30*30+1)=dRdl;
J(30*30+1,30*30+1)=-2*(lambda-lambda0);

% compute the dEds
b(30*30+1)=2*sqrt((lambda-lambda0)^2+norm(u-u0)^2);

% compute the duds and dlds
delta=-J\b;

%update the u and l
u=u+deltas*delta(1:30*30);
lambda=lambda+deltas*delta(30*30+1);

end


