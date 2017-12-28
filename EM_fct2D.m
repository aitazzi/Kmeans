function [mu,sigma,Pi,loglikelihood] = EM_fct2D(y,mu,sigma,Pi,NbIterations)

limite = 1;
N_em = length(Pi);
N = length(y);
p = zeros(N_em,N);

Pi_new = Pi;
sigma_carre_new = sigma;
mu_new = mu;
yp = y;
d=2;
while (limite < NbIterations)
    
    for n=1:N
        for k=1:N_em
            p(k,n) = Pi(k) * (2*pi)^(-d/2)/sqrt(det(sigma{k}))*exp(-0.5*(y(n,:)-mu(k,:))*inv(sigma{k})*(y(n,:)-mu(k,:))');
        end
        p(:,n) = p(:,n)/sum(p(1:N_em,n));
    end
    
    for k=1:N_em
        
        Pi_new(k) = sum(p(k,1:N))/N;
        for n=1:N
            yp(n,:) = y(n,:)*p(k,n);
        end
        mu_new(k,:) = sum(yp)/sum(p(k,1:N));
        
        for n=1:N
            zp(:,:,n) = p(k,n).*(y(n,:)-mu_new(k,:))'*(y(n,:)-mu_new(k,:));
            %*
        end
        sigma_carre_new{k} =  sum(zp,3)/ (sum(p(k,1:N)));
        
    end
    
    mu = mu_new;
    sigma = sigma_carre_new;
    Pi = Pi_new;
    
    limite = limite + 1;
end
likeli=zeros(n,N_em);
for i=1:N
    for k=1:N_em
        likeli(i,k)=p(k,i)*(log(Pi(k))-0.5*d*log(2*pi)-log(sqrt(det(sigma{k})))+(-0.5)*(y(i,:)-mu(k,:))*inv(sigma{k})*(y(i,:)-mu(k,:))');
    end
end
loglikelihood=sum(sum(likeli));
end