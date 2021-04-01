function[count]=solution(lambda_0)

%use the inital guess to generate solution iteratively

update = 0.1;

% parameter for initial guess
m = 2;
n=1;
A = -0.1;
u0 = initguessFunction(m,n,A);
[u0, J] = newtonMethod(u0, lambda_0);
lambda_1 = lambda_0 - update;
u1 = analyticContinuation(u0,-update,lambda_0, J);
u1 = newtonMethod(u1, lambda_1);
u_value = u1;
lambda = lambda_1;
count=zeros(2,10);
number = 1;
while lambda > 0 && norm(u_value) < 12
    [u_n, lambda_n, deltas] = alc(u_value, lambda, u0, lambda_0);
    [u_n, lambda_n] = augumentedNT(u_n, lambda_n, deltas, u0, lambda_0);
    u0 = u_value;
    u_value = u_n;
    lambda_0 = lambda;
    lambda = lambda_n;
    count(1, number) = norm(u_value);
    count(2,number) = lambda;
    number = number + 1;
    %draw the contour plot for the result
    if mod(number,20) == 0   
        dx=1/29;                     
        dy=1/29;                     
        [X,Y] = meshgrid(0:dx:1,0:dy:1);
        k=zeros(30);
        for i=1:30
        for j=1:30
            l=(j-1)*(30)+i;
            k(i,j)=u_value(l) ;
         end
         end
        Z = k;
        contour(X,Y,Z,'ShowText','on')
        colorbar
        xlabel('x')
        ylabel('y')
        title(['A=-0.1, 5*pi**2, m=2 n=1 ,lambda= ' num2str(lambda)])

        filename = ['b4 ',num2str(number),'.png'];
        saveas(gcf,filename)
        end
end

