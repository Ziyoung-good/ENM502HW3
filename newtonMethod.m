function [u, J] = newtonMethod(u, lambda)
%newton method code

% initialize the parameter
dx = 1/29;
dy = 1/29;
x = 0:dx:1;
y = 0:dy:1;
J = zeros(30*30);
R = zeros(30*30,1);
deltas = 1;

% %compute R and J until the result converge
while norm(deltas)>1e-3
for i=1:30
    for j=1:30
        index=(j-1)*(30)+i;
        if ((0<y(j))&&(y(j)<1)&&(0<x(i))&&(x(i)<1))
            R(index)=(u(index+1)-2*u(index)+u(index-1))/dx/dx+(u(index+30)-2*u(index)+u(index-30))/dy/dy+lambda*u(index)*(1+u(index));
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

%update the value
deltas=-J\R;
u=u+deltas;

end

end

