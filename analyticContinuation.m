function [u,lambda]=analyticContinuation(u,dl,lambda_0,J)
%ac function 


dx=1/29;                     
dy=1/29;                     
x=0:dx:1;                        
y=0:dy:1;
lambda=lambda_0+dl;
dRdl=zeros(30*30,1);

% compute the dRdlambda
for i=1:30
    for j=1:30
        index=(j-1)*(30)+i;
        if ((0<y(j))&&(y(j)<1)&&(0<x(i))&&(x(i)<1))
           dRdl(index)=u(index)*(1+u(index));
        else
           dRdl(index)=0;
         
        end
    end

end
% update the u value by using dl*dudl
u=u+(-J\dRdl)*dl;
end

