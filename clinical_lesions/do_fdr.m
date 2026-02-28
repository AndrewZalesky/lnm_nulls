function ind=do_fdr(p,Alpha)

[p_srt,ind_srt]=sort(p); 
ind=find(p_srt<Alpha*(1:length(p))/length(p));
if ~isempty(ind)
    ind=ind_srt(1:ind(end)); 
end
    