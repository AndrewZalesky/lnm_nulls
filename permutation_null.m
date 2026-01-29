function pval=permutation_null(sigma,N,mask,ind,Perms,alpha,rval,ind_inner,effect_size)

r_null=zeros(Perms,1);
for i=1:Perms
    img1_s=randn(N,N);
    img1_s=imgaussfilt(img1_s,sigma,'Padding','symmetric');
    img1_s(~mask)=0;

    img2_s=randn(N,N);
    img2_s=imgaussfilt(img2_s,sigma,'Padding','symmetric');
    img2_s(~mask)=0;

    img_tmp=randn(N,N);
    img_tmp=imgaussfilt(img_tmp,sigma,'Padding','symmetric');
    img_tmp(~mask)=0;
    %Preserve any increase in autocorrelation from the generative process
    img2_s=(1-alpha)*img2_s+alpha*img_tmp;

    if ~isempty(ind_inner)
        img_tmp=randn(N,N);
        img_tmp=imgaussfilt(img_tmp,sigma,'Padding','symmetric');
        img_tmp(~mask)=0;
        img1_s(ind_inner)=effect_size*img_tmp(ind_inner); 
        img2_s(ind_inner)=-effect_size*img_tmp(ind_inner); 
    end

    r_null(i)=corr(img1_s(ind),img2_s(ind));
end

pval=sum(r_null>=rval)/Perms; 