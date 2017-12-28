function [ cluster,center,distortion ] = k_means( data,K,init, affichage)

%% Main loop of the K-means algorithm (until convergence)
%affichage=1;
couleur = 'brmgky';
D = size(data,2);

x_min = double(min(data));  % to be sure tro work woth real numbers (not integers)
x_max = double(max(data));

nitermax = 200;   % maximum number of iterations
N = size(data,1) ;  % nb of pixels in the image
k=0;
cluster=zeros(N,1);
distance=zeros(N,1);
%Initialization
if init==0
    tmp = rand(1,K*D*2);  % to get K points in D dimensions
    center_new = repmat(x_min + (x_max-x_min),1,K).*tmp(1:(K*D)); % K centres
    center_new = reshape(center_new,K,D);
    center_old = repmat(x_min + (x_max-x_min),1,K).*tmp(K*D+1:end);
    center_old = reshape(center_old,K,D);
else
    center_old=data(randperm(N,K),:);
    center_new=data(randperm(N,K),:);
end

epsilon = 0.0000001*min(x_max-x_min); % stopping criterion on the evolution of centres
while k<nitermax && norm(center_old-center_new)/norm(center_old)>epsilon
    for i=1:N
        vector_distance=sqrt(sum(abs(repmat(data(i,:),K,1)-center_old).^2,2));
        [distMin,I]=min(vector_distance);
        cluster(i,:)=I;
        distance(i,:)=distMin^2;
    end
    center_old=center_new;
    for j=1:K
        Index=(cluster==j);
        center_new(j,:)=sum(data(Index,:),1)./sum(Index);  
    end  
    k=k+1;
end 
distortion=sum(distance);

center=center_new;
if affichage
    figure(1000)
    clf
    for j=1:K
        Index=(cluster==j);
        center_new(j,:)=sum(data(Index,:),1)./sum(Index);  
        set(gcf,'color','w')
        eval(['plot(data(Index,1),data(Index,2),''',couleur(j),'.'',''MarkerSize'',7)'])
        hold on
        eval(['plot(center_new(j,1),center_new(j,2),''',couleur(j),'p'',''MarkerSize'',24,''MarkerFaceColor'',''y'',''LineWidth'',2)'])
        grid on
    end
end
end
