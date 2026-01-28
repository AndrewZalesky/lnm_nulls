function [w,z]=dcsbm(GM)
%degree controlled stochatic block model (dcSBM)
%Generate a fully connected functional connectivity matrix with modules and
%log-normal degree distribution. The matrix is undirected, asymmetric and 
%not transitive. Connections weights are between 0 and 1, sampled from a 
%beta distribution. Future work: incorporate transitivity into the
%generative model. 
%w: connectivity matrix
%z: node allegiance 

%Probabilities
Bintra=GM.Bintra; %intra-edge 
Binter=GM.Binter; %inter-edge

%Number of modules
B=GM.B; 

%Number of nodes per module
M=GM.M;

sigma=GM.sigma; %lognormal distribution parameter
                %small sigma: homogeneous strength distribution
                %large sigma: strong hubs 

k=GM.k;         %dispersion in beta distribution for weights
                %large k: weights tightly coupled around expected value
                %small k, high variability

%Total number of nodes
N=B*M; 

%Node allegiance
z=zeros(N,1); 
for i=1:B
    z(1+(i-1)*M:i*M)=i; 
end

%Lognormal strength distribution
%sigma=0.2; %small sigma: homogeneous strength distribution
%           %large sigma: strong hubs 
theta=lognrnd(1,sigma,N,1); 

%Connection probabilities matrix 
if Bintra==Binter
    p=ones(N,N)*Bintra;
    p=p-diag(diag(p)); 
else
    p=zeros(N,N);
    for i=1:N
        for j=i+1:N
            if z(i)==z(j)
                p(i,j)=Bintra; p(j,i)=Bintra;
            else
                p(i,j)=Binter; p(j,i)=Binter;
            end
        end
    end
end

%Expected weight 
mu=theta*theta'.*p;
mu=mu/(max(mu(:))+eps);  

%k=100; %dispersion - large k, weights tightly coupled around mu
%       %           - small k, high variability

%Sample weights from beta distribution       
w=betarnd(mu*k,(1-mu)*k); 

%symmetrize
w=triu(w,1)+triu(w,1)';  