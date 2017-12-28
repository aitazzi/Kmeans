%% Load data
train=load('EMGaussian.data');
test=load('EMGaussian.test');
%% Plot training data
figure(1)
plot(train(:,1),train(:,2),'r*')
title('scatter of training data')
set(gcf,'color','w')
grid on
%% Plot testing data
figure(2)
plot(test(:,1),test(:,2),'b*')
title('scatter of training data')
set(gcf,'color','w')
grid on
%% Number of clusters
K=4;
%% Kmeans algorithm
display=1; % This parameter takes 1 if we want to display the kmeans result and 0 otherwise
init=1; %1 centres random among the data and 0 if random completely
[ cluster,center, distortion] = k_means(train,K,init,display);
%% Histogram of convergence ( nbr of clusters) and distortion
% Run if we want to compute the distrotions histogram
Nmax=1850;
cluster_num=zeros(Nmax,1);
%distortion_values=zeros(Nmax,1);
distortion_values=[];
display=0;
init=1;
for i=1:Nmax
    [ cluster,center,distortion ] = k_means(train,K,init,display);
    %cluster_num(i)=length(find(~isnan(center(:,1))));
    distortion_values=[distortion_values distortion];
end
figure(3)
set(gcf,'color','w')
hist(cluster_num,1:4)
title('histogram of the kmeans cluster number')
figure(4)
set(gcf,'color','w')
hist(distortion_values)
title('histogram of the values of distortions')
%% EM Algorithm Isotropic case
sigma_init=[1,1,1,1];
alpha_init=0.25.*ones(4,1);
NbIterations=10;
[mu,sigma,alpha,loglikelihood_iso] = EM_fct_isotropic(train,center,sigma_init,alpha_init,NbIterations);
Sigma={sigma(1)*eye(2),sigma(2)*eye(2),sigma(3)*eye(2),sigma(4)*eye(2)};
%% Plot the result
colors = {'r','b','g','k'};
markers = {'o','.','s','*'};
figure(3)
for k = 1:K
    idx = find(cluster == k);
    plot(train(idx,1),train(idx,2),sprintf('%s%s',colors{k},markers{k}));
    hold on
    plot_ellipse(Sigma{k},mu(k,:),colors{k});
    hold on 
    eval(['plot(mu(k,1),mu(k,2),''',colors{k},'p'',''MarkerSize'',24,''MarkerFaceColor'',''y'',''LineWidth'',2)'])
    set(gcf,'color','w')
end
title('EM isotropic case')
axis equal
%% EM general case
sigma_init={eye(2),eye(2),eye(2),eye(2)};
alpha_init=0.25.*ones(4,1);
NbIterations=100;
[mu,Sigma,alpha,loglikelihood] = EM_fct2D(train,center,sigma_init,alpha_init,NbIterations);
%% Plot the result
colors = {'r','b','g','k'};
markers = {'o','.','s','*'};
figure(3)
for k = 1:K
    idx = find(cluster == k);
    plot(train(idx,1),train(idx,2),sprintf('%s%s',colors{k},markers{k}));
    hold on
    plot_ellipse(Sigma{k},mu(k,:),colors{k});
    eval(['plot(mu(k,1),mu(k,2),''',colors{k},'p'',''MarkerSize'',24,''MarkerFaceColor'',''y'',''LineWidth'',2)'])
    set(gcf,'color','w')
end
title('EM general case')
axis equal
%% The loglikelihood for the testing data
% Isotropic
sigma_init=[1,1,1,1];
alpha_init=0.25.*ones(4,1);
NbIterations=100;
[mu_test,sigma_test,alpha_test,loglikelihood_iso_test] = EM_fct_isotropic(test,center,sigma_init,alpha_init,NbIterations);
%% 
%  general case
sigma_init={eye(2),eye(2),eye(2),eye(2)};
alpha_init=0.25.*ones(4,1);
NbIterations=100;
[mu_test,Sigma_test,alpha_test,loglikelihood_test] = EM_fct2D(test,center,sigma_init,alpha_init,NbIterations);

