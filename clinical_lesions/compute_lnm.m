function [lnm,lnm_dist,ii,jj]=compute_lnm(C,x,Trials)

%LNM
Nregs=sum(x,2); %number of separate regions representing each lesion
x_norm=x./repmat(Nregs,1,length(C)); 
lnm=x_norm*C; 
%lnm_map=mean(lnm); 

%LNM where only one region is sampled from each lesion
Nc=length(C);
lnm_dist=zeros(Trials,Nc); 
[ii,jj,~]=find(x); 
for n=1:Trials
    ind_region=accumarray(ii,jj,[],@(x) x(randi(length(x)))); % randomly sample one region per lesion
    lnm_dist(n,:)=mean(C(ind_region,:));
end
%lnm_dist_map=mean(lnm_dist); 