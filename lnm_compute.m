function [lnm,ind_sig,mod,ind_sig_ttest,null_dist]=lnm_compute(alpha,K,N,M,m,z,Assignment_Type,Perms,w,B,Null_Type)
%Assign lesions to regions based on alpha and Assignment_Type. Compute the
%LNM based on the connectivity matrix (w) and lesion locations. Perform
%statistical inference on the LNM using either the location or topology
%null model. Permutation testing is used with the maximum statistic in both
%cases, ensuring strong control of the family-wise error rate. More
%specifcially, the maximum value in the resampled LNM is stored for each
%permutation to build an empirical null distribution. Values for each
%node/region in the observed LNM are then compared to this null
%distribution, assuming a Type I error rate (corrected) of 5%. 
%lnm:     observed LNM, 1xN
%ind_sig: binary 1xN vector with a 1 positioned at each node/region
%         declared signifciant    
%mod:     binary Bx1 vector with a 1 position at each module declared
%         signficant under a Binomial test (alternative to LNM)
%
%Note that the Binomial test provides inference at the level of modules,
%whereas LNM enables inference at the higher resolution of nodes/regions. 

global w_null_topology

%Assign lesions
if Assignment_Type==2 || Assignment_Type==3
    %For observed data, types 2 and 3, both allow overlap in lesions 
    ind_lesion=lesion_assignment(alpha,K,N,M,m,z,2,[],[],[]); %Assignment type here must 
%                                                                     be either 1 or 2, never 3
elseif Assignment_Type==1
    ind_lesion=lesion_assignment(alpha,K,N,M,m,z,1,[],[],[]); 
end
%Compute LNM
lnm=ind_lesion'*w/K; 

s=sum(w)/N; %node strengths

%Neighbors for strength preservation
U=0.2; %top 20% of nodes matched in strength
un=ceil(U*(N-1));
neighbors = zeros(N, un);
for i = 1:N
    [~, ind_srt] = sort(abs(s(i) - s));
    neighbors(i,:) = ind_srt(2:un+1);
end

if strcmp(Null_Type,'Location')
    %Null model - lesion location randomization
    null_dist=zeros(Perms,1);
    for i=1:Perms
        ind_lesion_null=lesion_assignment(0,K,N,M,m,z,Assignment_Type,neighbors,ind_lesion,un);
        lnm_null=ind_lesion_null'*w/K;
        null_dist(i)=max(lnm_null); %max statistic, provides FWER correction
    end
elseif strcmp(Null_Type,'Topology')
    %Null model - generate network with no modular structure but similar degree distribution
    null_dist=zeros(Perms,1);
    for i=1:Perms
%         w_null=dcsbm(1,1,GM.B,GM.M,GM.sigma,GM.k); 
%         s_null=sum(w_null)/N; %node strengths;
%         w_null=w_null*mean(s)/mean(s_null); %re-scale to match mean degree
        lnm_null=ind_lesion'*w_null_topology{i}/K;
        null_dist(i)=max(lnm_null);
    end
end

null_dist_srt=sort(null_dist); 
critical_value=null_dist_srt(ceil(0.95*Perms)); 
ind_sig=lnm>critical_value; %signifciant nodes

%Binomial test (alternative to LNM in the special case where modules are
%known)
prop=1/B; %expected proportion of lesions per module under null of uniform assignment  
bonf=0.05/B; 
critical_value_bino=binoinv(1-bonf,K,prop)+1; 
mod=zeros(B,1); 
for i=1:B 
    if sum(ind_lesion(i==z))>critical_value_bino
        mod(i)=1; 
    end
end

%One-sample t-test on LNM, used for inference in some LNM work
%The null hupothesis is an LNM of zero
%Assume each lesion is a distinct observation
%This is a parametric test. Non-parametric sign flipping could also be used
cnt=1; 
subsamp=zeros(K,N); 
for i=1:N
    if ind_lesion(i)>0
        for j=1:ind_lesion(i)
            subsamp(cnt,:)=w(i,:);
            cnt=cnt+1; 
        end
    end
end
[~,pval]=ttest(subsamp);
ind_sig_ttest=pval<0.05/N; %Bonf correction 
