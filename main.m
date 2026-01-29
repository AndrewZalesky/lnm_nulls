% Evaluate null models for a bare-bones LNM. Functional connectivity
% matrices generated from a degree-corrected stochastic block model
% (dcSBM). Two null models are considered: (i) lesion location
% randomization and (ii) network topology rewiring/regeneration. The ground
% truth LNM is defined by a single target network module. A proportion of
% lesions (e.g., alpha=0.8) are distributed across nodes comprising the
% target module. Remaining lesions are distributed uniformly across the
% other modules. When alpha=0, all lesions are distributed uniformly,
% irrespective of modular structure, enabling evaluation of false positive
% rate control of the LNM method. When alpha>0, the statistical power to
% detect an effect under the two null models is evaluated. For each LNM,
% the maximum statistic is used to construct the null distribution across
% thousands of lesion randomizations or network rewirings/regenerations,
% ensuring strong control of the family-wise error rate (FWER). Note that
% all evaluations are limited to LNMs defined by a single predefined target
% module. In this simple case where modular composition is known, as an
% alternative to the LNM method, we can also simply test whether lesions
% are more likely to reside in a given module than expected under a uniform
% distribution of lesions, using a Binomial distribution. However, in
% practice, modular organization may not be uniquely definable (i.e.,
% modules can be defined across different scales, hierarchically,
% overlapping, etc.). An advantage of the LNM method  is that it does not
% make any assumptions about the topology underlying the LNM.
%Note that all simulations are based on a regional LNM model. Extending the
%topological null to voxel-resolution LNM may be computationally
%challenging.
clear all
close all

%Add colors
addpath('./cbrewer/'); 
cmap=cbrewer('qual','Set3',12);
cmap_matrix=cbrewer('seq','Reds',256); 

%Parameters for dcSBM
Bintra=0.6; %probability of intra-modular connection
Binter=0.3; %probability of inter-modular connection
B=10;       %number of modules
M=20;       %number of nodes per module
sigma=0.2;  %lognormal distribution parameter
            %small sigma: homogeneous strength distribution
            %large sigma: strong hubs 
k=10;       %dispersion in beta distribution for weights
            %large k: weights tightly coupled around expected value
            %small k, high variability
N=B*M;      %number of nodes
GM=struct('Bintra',Bintra,...
          'Binter',Binter,...
          'B',B,...
          'M',M,...
          'sigma',sigma, ...
          'k',k,...
          'N',N); 

%LNM parameters
K=100;      %number of lesions 
alpha=0.8;  %proportion of lesions residing within selected module
            %remianing lesions distributed uniformly across other modules
            %this is an effect size
m=5;        %select a module (arbitrary)
%Assignment of lesions to regions
%1 or 2
%1: region/voxel assigned at most one lesion
%2: region/voxel can be assigned multiple lesions
Assignment_Type=2;
LNM=struct('K',K,...
           'alpha',alpha,...
           'm',m,...,
           'Assignment_Type',Assignment_Type); 

%Statistical parameters
Perms=5000; %number of permutations for null testing
%Type of null model
Null_Type='Location'; %randomise lesion locations
%Null_Type='Topology'; %randomise (re-generate) topology without modules
%Number of trials for averaging
Trials=1000;  %set to 10000 for generating figures 

%Generate NxN symmetric functional connectivity matrix
%Nodes could represent voxels or regions 
[w,z]=dcsbm(GM); 
%w: connectivity matrix
%z: Nx1 vector of node to module assignments, modules labeled 1,...,B

s=sum(w)/N; %node strengths

fprintf('Null type: %s\n',Null_Type); 
fprintf('Defaults: alpha=%0.1f, K=%d, Trials=%d, Perms=%d\n',alpha,K,Trials,Perms);


%Precompute topology rewiring/regeneration to avoid rewiring/regeneration
%for each trial
global w_null_topology %global variable to store randomized topologies
w_null_topology=cell(Perms,1); 
precompute_topology_randomization(s,GM,Null_Type,Perms); 
%Note that the node strength distribution in the null samples is matched to 
%the single network sample generated above. 
%Each trial below will generate a new sample with a
%potentially slightly different strength distribution. To ensure perfect
%matching of strength distributions, no precomputation should be perfomed.
%But this slows exceution considerably. 

%%
%Figure: Functional connectivity matrix and node strength distribution
hf1=figure; hf1.Position=[100,100,400,400]; hf1.Color='w';
               imagesc(w,[0,1]); colormap(cmap_matrix); axis equal; axis off; %title('Connectivity matrix');
               %colorbar;
                ha11=gca; ha11.FontSize=20; 
hf1a=figure; hf1a.Position=[100,100,300,150]; hf1a.Color='w';  
                [f,xi]=ksdensity(s,'Bandwidth',0.2,'Support','Positive');
                plot(xi,f,'Linewidth',2); %title('Node strength distribution'); 
                ha12=gca; ha12.FontSize=20; ha12.LineWidth=3;
%%
if strcmp(Null_Type,'Topology')
    %%
    %Figure: Functional connectivity matrix under null hypothesis and node strength distribution
hf1null=figure; hf1null.Position=[100,100,400,400]; hf1null.Color='w';
               imagesc(w_null_topology{1},[0,1]); colormap(cmap_matrix); axis equal; axis off; %title('Connectivity matrix');
               %colorbar;
                ha11null=gca; ha11null.FontSize=20; 
hf1anull=figure; hf1anull.Position=[100,100,300,150]; hf1anull.Color='w';  
                s_null=sum(w_null_topology{1})/N; 
                [f,xi]=ksdensity(s_null,'Bandwidth',0.2,'Support','Positive');
                plot(xi,f,'Linewidth',2); %title('Node strength distribution'); 
                ha12null=gca; ha12null.FontSize=20; ha12null.LineWidth=3;
    %%
end

%Evaluate as a function of alpha, fixed K
Mode='test_alpha'; 
%Evaluate as a function of K, fixed alpha
%Mode='test_K';
%Evaluate a single value of alpha and K, and
%generate a detailed figure 
%Mode='test_single_case'; 

if strcmp(Mode,'test_alpha')
    fprintf('Evaluating as a function of alpha...\n'); 
    %Evaluate as a function of alpha
    alpha_rng=[0:0.05:1]; J=length(alpha_rng); 
    %Rates of interest - see below for description
    fwer=zeros(J,1); 
    power=zeros(J,1); 
    strong_fwer=zeros(J,1); 
    tpr_binomial=zeros(J,1);
    fpr_module_binomial=zeros(J,1);
    tpr_ttest=zeros(J,1); 
    frst=0; 
    for i=1:J
        alpha_val=alpha_rng(i);
        for n=1:Trials
            [w,z]=dcsbm(GM); 
            s=sum(w)/N; %node strengths
            [lnm,ind_sig,mod,ind_sig_ttest,null_dist]=lnm_compute(alpha_val,K,N,M,m,z,Assignment_Type,Perms,w,B,Null_Type);
            %LNM
            fwer(i)=fwer(i)+any(ind_sig); %fwer but only when alpha=0 (null condition)
            power(i)=power(i)+any(ind_sig(z==m)); %power when alpha>0
            strong_fwer(i)=strong_fwer(i)+any(ind_sig(z~=m)); %strong fwer under the alternative hypothesis                                        
            %Binomial test
            %Prob of detecting an effect in at least one module  
            tpr_binomial(i)=tpr_binomial(i)+any(mod);
            %Prob of one or more modules detected other than the selected
            %target module
            fpr_module_binomial(i)=fpr_module_binomial(i)+any(mod([1:m-1,m+1:B])); 
            %Note that resolution of inference is different between LNM and
            %Binomial. Former is at nodes, latter is at modules. 
            tpr_ttest(i)=tpr_ttest(i)+any(ind_sig_ttest); 
        end
        show_progress(i*2,J*2,frst); frst=1; 
    end
    fwer=fwer/Trials; power=power/Trials; strong_fwer=strong_fwer/Trials;  
    tpr_binomial=tpr_binomial/Trials; fpr_module_binomial=fpr_module_binomial/Trials; 
    tpr_ttest=tpr_ttest/Trials; 
    %%
    %Figure: Evaluation of alpha
    hfA=figure; hfA.Position=[100,100,1000,400]; hfA.Color='w';
    subplot(1,2,1); plot(alpha_rng,[power,fwer]); legend('Power','FWER'); xlabel('alpha');
    subplot(1,2,2); plot(alpha_rng,strong_fwer); legend('FWER Strong'); xlabel('alpha'); 
    %%
    save(strcat(['alpha_data_',Null_Type,'_Kval',num2str(K),'.mat']),...
                                                  'alpha_rng','fwer','power','strong_fwer','tpr_binomial',...
                                                  'fpr_module_binomial','tpr_ttest');

elseif strcmp(Mode,'test_K')
    fprintf('Evaluating as a function of K...\n'); 
    %Evaluate as a function of K
    K_rng=[1,10,20,50,100,120]; J=length(K_rng); 
    fwer=zeros(J,1); 
    power=zeros(J,1); 
    strong_fwer=zeros(J,1); 
    tpr_binomial=zeros(J,1);
    fpr_module_binomial=zeros(J,1);
    tpr_ttest=zeros(J,1); 
    frst=0; 
    for i=1:J
        K_val=K_rng(i);
        for n=1:Trials
            [w,z]=dcsbm(GM); 
            s=sum(w)/N; %node strengths
            [lnm,ind_sig,mod,ind_sig_ttest,null_dist]=lnm_compute(alpha,K_val,N,M,m,z,Assignment_Type,Perms,w,B,Null_Type);
            fwer(i)=fwer(i)+any(ind_sig); %fwer but only when alpha=0 (null condition)
            power(i)=power(i)+any(ind_sig(z==m)); %power when alpha>0
            strong_fwer(i)=strong_fwer(i)+any(ind_sig(z~=m)); %strong fwer under the alternative hypothesis                                        
            %Binomial test
            %Prob of detecting an effect in at least one module  
            tpr_binomial(i)=tpr_binomial(i)+any(mod);
            %Prob of one or more modules detected other than the selected
            %target module
            fpr_module_binomial(i)=fpr_module_binomial(i)+any(mod([1:m-1,m+1:B])); 
            %Note that resolution of inference is different between LNM and
            %Binomial. Former is at nodes, latter is at modules. 
            tpr_ttest(i)=tpr_ttest(i)+any(ind_sig_ttest);  
        end
        show_progress(2*i,2*J,frst); frst=1; 
    end
    fwer=fwer/Trials; power=power/Trials; strong_fwer=strong_fwer/Trials;  
    tpr_binomial=tpr_binomial/Trials; fpr_module_binomial=fpr_module_binomial/Trials; 
    tpr_ttest=tpr_ttest/Trials;
    %%
    %Figure: Evaluation of K
    hfA=figure; hfA.Position=[100,100,1000,400]; hfA.Color='w';
    subplot(1,2,1); plot(K_rng,[power,fwer]); legend('Power','FWER'); xlabel('Number of lesions');
    subplot(1,2,2); plot(K_rng,strong_fwer); legend('Strong FWER');
                    xlabel('Number of lesions'); 
    %%
    save(strcat(['K_data_',Null_Type,'_Kval',num2str(K),'.mat']),...
                                              'K_rng','fwer','power','strong_fwer','tpr_binomial',...
                                              'fpr_module_binomial','tpr_ttest'); 

elseif strcmp(Mode,'test_single_case')
    %Generate detailed figure for default parameters, specfied above 

    [lnm,ind_sig,mod,~,null_dist]=lnm_compute(alpha,K,N,M,m,z,Assignment_Type,Perms,w,B,Null_Type);

    %%
    %Figure 2
    hf2=figure; hf2.Position=[150,150,700,300]; hf2.Color='w';
    stem(s,'.k','LineWidth',1.5,'MarkerSize',10); hold;
    stem(find(~ind_sig),lnm(~ind_sig),'.r','LineWidth',1.3,'MarkerSize',10);
    stem(find(ind_sig),lnm(ind_sig),'*r','LineWidth',1.3,'MarkerSize',10);
    gap=0.05;
    min_val=min([lnm,s])-gap; max_val=max([lnm,s]);
    ylim([min_val,0.6]);
    ha2=gca; ha2.XTick=[]; ha2.LineWidth=3;  ha2.Box='off'; 
    %ylabel('Connectivity');
    for i=1:B
        x1=(i-1)*M+1-0.5;
        x2=i*M+0.5;
        y1=min_val;
        y2=min_val+gap;
        x=[x1 x1 x2 x2];
        y=[y1 y2 y2 y1];
        p=patch(ha2,x,y,cmap(i,:));
        p.FaceAlpha=0.8;
        p.EdgeColor='none';
    end
    hl=legend('Node strength','LNM');
    hl.EdgeColor=[1,1,1];
    hl.Location='northwest'; 
    xlim([1,N]);
    ha2.FontSize=20;

%Null distribution histogram 
%     hfnull=figure; hfnull.Position=[100,100,250,300]; hfnull.Color='w';  
%         [f,xi]=ksdensity(null_dist,'Bandwidth',0.2,'Support','Positive');
%          plot(xi,f,'Linewidth',2); %title('Node strength distribution'); 
%          h=histogram(null_dist,40,'FaceColor',cmap(1,:)); 
%          hanull=gca; hanull.FontSize=20; hanull.LineWidth=3;
%          hanull.XTick=[0.3,0.4,0.5,0.6]; 
%             xticks([0.3,0.35,0.4,0.45]);
    %%
end