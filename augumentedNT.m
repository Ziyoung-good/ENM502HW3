function [u,lambda]=augumentedNT(u,lambda,ds,u0,lambda0)
%Full Newton’s Method with an Augmented Residual/Jacobian
                                 
dx=1/29;                     
dy=1/29;                     
x=0:dx:1;                        
y=0:dy:1; 
R=zeros(30*30+1,1);
delta=1;
J=zeros(30*30+1);
dRdl=zeros(30*30,1);
dEdu=zeros(30*30,1);

%compute R and J which is similar to the newton method
while norm(delta)>1e-3
for i=1:30
    for j=1:30
        index=(j-1)*(30)+i;
        dEdu(index)=-2*(u(index)-u0(index));
        if ((0<y(j))&&(y(j)<1)&&(0<x(i))&&(x(i)<1))
            R(index)=(u(index+1)-2*u(index)+u(index-1))/dx/dx+(u(index+30)-2*u(index)+u(index-30))/dy/dy+lambda*u(index)*(1+u(index));
            dRdl(index)=u(index)*(1+u(index));
            J(index,index)=-2/(dx*dx)-2/(dy*dy)+lambda*(1+2*u(index));
            J(index,index-1)=1/dx/dx;
            J(index,index+1)=1/dx/dx;
            J(index,index+30)=1/dy/dy;
            J(index,index-30)=1/dy/dy;
        else
           
            R(index)=u(index);
            J(index,index)=1;
        end
    end
end

%compute the J
J(30*30+1,1:30*30)=dEdu;
J(1:30*30,30*30+1)=dRdl;
J(30*30+1,30*30+1)=-2*(lambda-lambda0);

%compute the delta to udpate the u and lambda
R(30*30+1)=ds^2-norm(u-u0)^2-(lambda-lambda0)^2;


%update the u and lambda
delta=-J\R;
u=u+delta(1:30*30);
lambda=lambda+delta(30*30+1);


end
end

