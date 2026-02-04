function ind_lesion=lesion_assignment(alpha,K,N,M,m,z,Assignment_Type,neighbors,ind_lesion_observed,un)
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
%           3: Regions can be assigned to multiple lesions and matching for
%           randomization approximately preserves node strength
%s:         node strength observed (only relevant for spatial null)
%ind_lesion_observed: lesions in observed data 
%Note that s and ind_lesion_observed are only used for Assignment_Type==3. 

k_module=floor(alpha*K); 
k_other=K-k_module; 

ind_lesion=zeros(N,1);
if alpha==0
    %Assign lesions randomly, blind to modular structure; 
    if Assignment_Type==1
        ind_lesion(randperm(N,K))=1;
    elseif Assignment_Type==2
        idx=randi(N,K,1);
        ind_lesion=ind_lesion + accumarray(idx(:),1,[N,1]);
    elseif Assignment_Type==3
        %assign in a way that approximately matches node strength
        %Nodes are sorted base on how well they are matched in strength of
        %the target node. The top U% of nodes are then candidates for
        %repositioning. If U is small, strength matching is strong, but
        %diversity in permutations is low. If U is high, matching in
        %strength is more approximate, but more diversity in randomization.

        %%%Too slow 
%         for i=1:N
%             for j=1:ind_lesion_observed(i)
%                 [~,ind_srt]=sort(abs(s(i)-s)); 
%                 ind_srt=ind_srt(2:un+1); %remove the node under consideration
%                                         %to avoid repositing to same node
%                 tmp=ind_srt(randperm(un,1));
%                 ind_lesion(tmp)=ind_lesion(tmp)+1; 
%                 %this allows the same region to be assigned multiple lesion
%                 %as per Assignment_Type 2
%             end
%         end
        %%%%
        %Faster version
        for i = 1:N
            k = ind_lesion_observed(i);
            if k > 0
                picks = neighbors(i, randi(un, k, 1));
                ind_lesion = ind_lesion + accumarray(picks', 1, [N,1]);
            end
        end
        %%%%
        %close all; 
        %figure; stem(ind_lesion_observed);
        %figure; stem(ind_lesion);
        %pause


    end
else
    %Assign lesions biased towards module m
    if Assignment_Type==1 
        %Option 1: regions can be assigned at most one lesion
        ind=find(m==z); 
        ind_lesion(ind(randperm(M,k_module)))=1;
        ind=find(m~=z);
        ind_lesion(ind(randperm(N-M,k_other)))=1;
    elseif Assignment_Type==2 
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
