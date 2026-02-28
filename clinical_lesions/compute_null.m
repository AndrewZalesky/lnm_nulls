function [lnm_dist_null,d]=compute_null(ii,jj,C,s,ind_srt,Nc,Perms,Nprop)

%Null - approximate strength matching 
frst=0; 
d=zeros(Perms,1);
lnm_dist_null=zeros(Perms,Nc); 
for n=1:Perms
    ind_region=accumarray(ii,jj,[],@(x) x(randi(length(x)))); % randomly sample one region per lesion, exactly as per real LNM
    %Note that the sample drawn above may be different from the sample
    %drawn in the observed LNM
    ind_tmp=sub2ind(size(ind_srt),ind_region,randi([1,Nprop],length(ind_region),1)); 
    ind_region_null=ind_srt(ind_tmp); %sample a region from each row approximately matched in strength
    d(n)=mean(abs(s(ind_region)-s(ind_region_null))); %mean distance between observed strength and null sample strength 
    lnm_dist_null(n,:)=mean(C(ind_region_null,:));
    show_progress(n,Perms,frst);frst=1; 
end