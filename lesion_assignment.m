function ind_lesion=lesion_assignment(alpha,K,N,M,m,z,Assignment_Type)
%Assign lesions to regions/nodes
%alpha:     Proportion of lesions assigned to a given module
%           When alpha=0, lesions assigned randomly
%K:         Number of lesions
%N:         Number of nodes
%M:         Number of modules
%m:         Selected module, ignored when alpha=0
%z          Module index, ignored when alpha=0
%Assignment_Type:
%           1: Regions can be assigned at most one lesion
%           2: Regions can be assigned multiple lesions 

k_module=floor(alpha*K); 
k_other=K-k_module; 

ind_lesion=zeros(N,1);
if alpha==0
    %Assign lesions randomly, blind to modular structure; 
    if Assignment_Type==1
        ind_lesion(randperm(N,K))=1;
    else
        idx=randi(N,K,1);
        ind_lesion=ind_lesion + accumarray(idx(:),1,[N,1]);
    end
else
    %Assign lesions biased towards module m
    if Assignment_Type==1 
        %Option 1: regions can be assigned at most one lesion
        ind=find(m==z); 
        ind_lesion(ind(randperm(M,k_module)))=1;
        ind=find(m~=z);
        ind_lesion(ind(randperm(N-M,k_other)))=1;
    else 
        %Option 2: regions can be assigned multiple lesions
        ind=find(m==z); 
        tmp=randi(M,k_module,1);
        idx=ind(tmp);
        ind_lesion=ind_lesion + accumarray(idx(:),1,[N,1]);
      
        %ind=find(m~=z); %Don't allow assignment to target module 
        %tmp=randi(N-M,k_other,1);

        ind=1:N;  %Allow assignment to all modules including target
        tmp=randi(N,k_other,1);

        idx=ind(tmp);
        ind_lesion=ind_lesion + accumarray(idx(:),1,[N,1]);
    end
end
