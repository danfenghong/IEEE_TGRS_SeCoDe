function [A, X] = SPCLSU_ADMM(Y,A, maxiter)

epsilon = 1e-6;
iter=0;

[l,N]=size(Y);
[l,k]=size(A);

P=zeros(l,k);
Q=zeros(k,N);

lamda1=zeros(size(P));
lamda2=zeros(size(Q));

stop = false;
mu=1e-3;
rho=1.5;
mu_bar=1E+6;

while ~stop && iter < maxiter+1
    
    iter=iter+1;
    
    X=(A'*A + mu*eye(size(A'*A)))\(A'*Y +mu * Q+ lamda2); 

    A = (Y * X' + mu * P +lamda1)/(X*X'+mu * eye(size(X*X')));
        
    P=max(A-lamda1/mu,0);  
    Q=max(X-lamda2/mu,0); 


   lamda1=lamda1+mu*(P-A);
   lamda2=lamda2+mu*(Q-X);

   mu=min(mu*rho,mu_bar);
   
    r_P=norm(P-A,'fro');
    r_Q=norm(Q-X,'fro');

    if r_P<epsilon&&r_Q<epsilon
            stop = true;
            break;
    end
end

end