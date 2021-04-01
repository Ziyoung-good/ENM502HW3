function [u] = initguessFunction(m, n, A)
%Get the initial guess function for u which correspond to the step one in
%problem setup part
% 30*30 grid
 
dx=1/29;                     
dy=1/29;                     
x=0:dx:1;                        
y=0:dy:1; 

u=zeros(30*30,1);
for i=1:30
    for j=1:30
        point=(j-1)*(30)+i;
        u(point) = A*sin(m*pi*x(i))*sin(n*pi*y(j));
            
    end

end

