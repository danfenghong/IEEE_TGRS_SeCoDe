function [A, X] = NMF_ADMM(Y, Z,alfa, beta, maxiter)

epsilon = 1e-6;
iter=0;

[l,N]=size(Y);
[k,N]=size(Z);

P=zeros(l,k);
Q=zeros(size(Z));
G=zeros(size(Z));

lamda1=zeros(size(P));
lamda2=zeros(size(Q));
lamda3=zeros(size(Z));

X = Z;

stop = false;
mu=1e-3;
rho=1.5;
mu_bar=1E+6;

while ~stop && iter < maxiter+1
    
    iter=iter+1;
    
    A = (Y * X' + mu * P +lamda1)/(X*X'+mu * eye(size(X*X')));
      
    X=(A'*A + (2*mu+beta)*eye(size(A'*A)))\(A'*Y + beta * Z +mu * Q+ lamda2 + mu * G + lamda3); 
    
    G=max(abs(X-lamda3/mu)-(alfa/mu),0).*sign(X-lamda3/mu); 
        
    P=max(A-lamda1/mu,0);  
    Q=max(X-lamda2/mu,0); 


   lamda1=lamda1+mu*(P-A);
   lamda2=lamda2+mu*(Q-X);
   lamda3=lamda3+mu*(G-X);

   mu=min(mu*rho,mu_bar);
   
    r_P=norm(P-A,'fro');
    r_Q=norm(Q-X,'fro');
    r_G=norm(G-X,'fro');

    if r_P<epsilon&&r_Q<epsilon&&r_G<epsilon
            stop = true;
            break;
    end
end

% X = Q;
% A = P;
end