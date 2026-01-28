function precompute_topology_randomization(s,GM,Null_Type,Perms)

global w_null_topology

%Set inter and intra probabilities to be equal, thus avoding any modular
%structure and generating a non-modular network with log-normal degree
%distribution.
GM.Binter=1; 
GM.Bintra=1; 

N=GM.N; 

if strcmp(Null_Type,'Topology')
    fprintf('Precomputing null for topology...\n');
    w_null_topology=cell(Perms,1); frst=0;
    for i=1:Perms
        w_null=dcsbm(GM);
        s_null=sum(w_null)/N; %node strengths;
        w_null=w_null*mean(s)/mean(s_null); %re-scale to match mean degree
        w_null_topology{i}=w_null;
        show_progress(i,Perms,frst); frst=1;
    end
else 
    w_null_topology=[]; 
end