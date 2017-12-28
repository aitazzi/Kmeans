function [mu,sigma,Pi,loglikelihood] = EM_fct_isotropic(y,mu,sigma,Pi,NbIterations)
limite = 1;
N_em = length(Pi);
N = length(y);
p = zeros(N_em,N);

Pi_new = Pi;
sigma_carre_new = sigma;
mu_new = mu;
yp = y;
d=2;
while  (limite < NbIterations) 
    
    % Computing tau
    for n=1:N
        for k=1:N_em
            p(k,n) = Pi(k) * (2*pi)^(-d/2)/(sigma(k))*exp((-0.5/sigma(k))*(y(n,:)-mu(k,:))*(y(n,:)-mu(k,:))');
        end
        p(:,n) = p(:,n)/sum(p(1:N_em,n));
    end
    
    % m STEP
    for k=1:N_em
        % Estimation of the Pi
        Pi_new(k) = sum(p(k,1:N))/N;
        % Estimation of the centers mu
        for n=1:N
            yp(n,:) = y(n,:)*p(k,n);
        end
        mu_new(k,:) = sum(yp)/sum(p(k,1:N));
        % Estimation of the covariance matrices sigma
        for n=1:N
            zp(n) = p(k,n).*(y(n,:)-mu_new(k,:))*(y(n,:)-mu_new(k,:))';
        end
        sigma_carre_new(k) =  sum(zp)/ (2*sum(p(k,1:N)));
    end
    
    mu = mu_new;
    sigma = sigma_carre_new;
    Pi = Pi_new;
    limite = limite + 1;
end
likeli=zeros(n,N_em);
for i=1:N
    for k=1:N_em
        likeli(i,k)=p(k,i)*(log(Pi(k))-0.5*d*log(2*pi)-log(sigma(k))+(-0.5/sigma(k))*(y(i,:)-mu(k,:))*(y(i,:)-mu(k,:))');
    end
end
loglikelihood=sum(sum(likeli));
end