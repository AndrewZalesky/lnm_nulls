clear all
close all

%Number of times regions are sampled within a lesion
Trials=500;

%Number of permutations for null distribution
Perms=5000;

Prop=0.2; %only allow permutations within the top Prop proportion of nodes 
          %closest in strength to the actual node

%p-value tail type {'positive','negative','twotail'}
%test='positive';
%test='negative';
test='twotail'; 

%Lesion sets to consider
L_set=[1,2]; 
%L_set=1;
%L{1}='emo'; 
%Lesion sets defined in Van den Heuvel et al                 
L{1}='SyntheticLesions';   %ns (expected)
L{2}='ADDICTION';          %ns             
L{3}='AGENCY';             %ns
L{4}='APHASIA_recovery';   %reject null    
L{5}='EPILEPSY';           %reject null                        
L{6}='LESYMAP';            %reject null    
L{7}='MIGRAINE';           %very few lesions, exclude           
L{8}='PARKINSONISM';       %ns
L{9}='NEGLECT';            %ns            
L{10}='OCD';               %ns
L{11}='STUTTERING';        %reject null                        

%Load lesions maps used by van den Heuvel et al 
%L is a list of lesions sets for each disorder
%all_lesions_atlas{i} is a lesions x regions binary matrix indicating
%affected regions for the disorder L{i} 
%all_lesion_atlas was mapped using map_voxel2schaeferMelbourne.m, as
%implemented in Example1_LNM_step123.m  
%No thresholding 
load('./data/lesion_maps.mat');
%load('./data/emo_lesions_map.mat');

%Load list of region indexes and names
load('./data/ctab.mat');

%Load GSP1000 regional connectivity matrix
%C is the N x N connectivity matrix used by van den Heuvel et al 
%Croi is a list or region names
load('./data/connectivity_matrix.mat'); 
%C=abs(C); 

%Read atlas volume 
[~,atlas_vol]=read('schaefer1000-yeo7+melbourne_MNI152_2mm.nii'); 
[~,atlas_vol_dil]=read('schaefer1000-yeo7+melbourne_MNI152_2mm_dilated.nii');

%%%MAPPINGS BETWEEN ATLAS REGIONS AND CONNECTIVITY MATRIX ROWS
%Establish mappings from atlas regions to connecitvity matrix and converse
atlas_vals=setdiff(unique(atlas_vol(:)),0); %unique region indexes
atlas2C=zeros(length(atlas_vals),1); 
%atlas2C(i) is the index of the ith atlas region in the connectivity matrix
for i=1:length(atlas_vals)
    ind_tmp=find(atlas_vals(i)==table2array(ctab(:,1)));
    atlas2C(i)=find(strcmp(Croi,ctab{ind_tmp,2})); 
end
C2atlas=zeros(length(atlas_vals),1); 
%C2atlas(i) is the index of the ith column/row of the connectivity matrix
%in the list of atlas regions (i.e., atlas_vals)
for i=1:length(atlas_vals)
    C2atlas(i)=find(i==atlas2C); 
end

%Number of regions
Nc=length(atlas_vals); 

%%%REMOVE DISTANCE EFFECTS FROM FUNCTIONAL CONNECTIVITY
%Load functional connectivity and regress out distance effects from 
%functional connectivity 
if ~exist('./data/connectivity_matrix_corrected.mat','file')
    fprintf('Residualizing distance from FC...\n'); 
    coor=zeros(length(atlas_vals),3);
    for i=1:length(atlas_vals)
        [ii,jj,kk]=find(atlas_vals(i)==atlas_vol_dil);
        coor(C2atlas(i),:)=[mean(ii),mean(jj),mean(kk)];
    end
    %Distance between all paris of regions
    D=squareform(pdist(coor)); %Euclidean is default 

    ind_upper=find(triu(ones(Nc,Nc),1));
    M=Nc*(Nc-1)/2;
    design=zeros(M,6); %design matrix for regression
    design(:,1)=1; %intercept
    design(:,2)=D(ind_upper); %distance effect
    design(:,3)=C(ind_upper)<0; %intercept for -ve FC
    design(:,4)=(C(ind_upper)<0).*D(ind_upper); %distance effect for -ve FC
    design(:,5)=D(ind_upper).^2; %distance effect squared
    design(:,6)=(C(ind_upper)<0).*(D(ind_upper).^2); %squared effect for -ve FC 
    %Regress distance from FC 
    [b,dev,stats]=glmfit(design,C(ind_upper),'normal','Constant','off');
    Cd=zeros(Nc,Nc);
    Cd(ind_upper)=stats.resid;
    Cd=Cd+Cd';
    hf=figure; hf.Position=[100,100,800,500]; hf.Color='w';
        subplot(2,2,1); dscatter(D(ind_upper),C(ind_upper)); xlabel('Distance'); ylabel('FC');
        subplot(2,2,2); dscatter(D(ind_upper),Cd(ind_upper)); xlabel('Distance'); ylabel('Distance Residualized');
        subplot(2,2,3); imagesc(C); title('FC'); 
        subplot(2,2,4); imagesc(Cd); title('Distance residualized FC'); 
    C=Cd; %replace FC with residualized FC
    save('./data/connectivity_matrix_corrected.mat','C','Croi')
else
    fprintf('Loading distance residualized FC...\n'); 
    load('./data/connectivity_matrix_corrected.mat')
end

%%%REGIONAL STRENGTH
s=sum(C)/length(C); %regional strength 
Nc=length(C);       %number of regions


%%%STRENGTH MATCHING MATRIX FOR NULL MODEL
%Strength matching matrix 
%ind_srt(i,:) is a an ordered list of regions indexes, ind_srt(i,1) is the region that
%is closest in strength to node i, ind_srt(i,2) is the region that is 2nd
%closest in strength to node i, etc. 
%To do: nodes closest in strength may also be neighbors, potentially add a
%spatial constraint. 
ind_s=zeros(Nc,Nc); 
d=abs(repmat(s,Nc,1)'-repmat(s,Nc,1));
[~,ind_srt]=sort(d,2);
ind_srt=ind_srt(:,2:ceil(Prop*Nc)); %start from 2 to exclude the reference node 
Nprop=size(ind_srt,2); %number of candiates

r=zeros(length(L_set),1);
ind=cell(length(L_set),1);
p=cell(length(L_set),1);
lnm_store=zeros(length(L_set),Nc); 
lnm_map=cell(length(L_set),1);
x_les=cell(length(L_set),Nc); 
lnm_full=cell(length(L_set),Nc);
for n=1:length(L_set)
    fprintf('Processing %s...\n',L{L_set(n)}); 
    x=all_lesions_atlas{L_set(n)}; 
    %Remove rows with no lesions
    x=x(sum(x,2)~=0,:); 
    x_les{n}=x; 
    fprintf('Number of lesions: %d\n',size(x,1));
    
    %Compute LNM
    %lnm is a lesions x regions matrix (classical LNM)
    %lnm_dist is a Trials x regions matrix
    %Each row is a connecitvity profile for a lesion (lnm) or single region
    %within a lesion (lnm_dist)
    [lnm,lnm_dist,ii,jj]=compute_lnm(C,x,Trials);
    %Correlation between classical and distributional LNM
    r(n)=corr(mean(lnm)',mean(lnm_dist)'); 
    lnm_store(n,:)=mean(lnm_dist); 
    lnm_full{n}=lnm; 

    %Null model
    %lnm_dist_null is Perms x regions matrix
    %d is Perms x 1 and used to assess accuracy of strength matching
    [lnm_dist_null,d]=compute_null(ii,jj,C,s,ind_srt,Nc,Perms,Nprop);

    

    %Compute (uncorrected) p-value for each region 
    %This invovles a double integral over samples of observed statistic and null 
    p{n}=zeros(Nc,1); %uncorrected p
    for i=1:Nc
        if strcmp(test,'twotail')
            p{n}(i)=...
                    (1+sum(sum(repmat(abs(lnm_dist(:,i)),1,Perms)...
                    <=repmat(abs(lnm_dist_null(:,i)),1,Trials)')))/(Trials*Perms+1); 
        elseif strcmp(test,'positive')
            p{n}(i)=...
                    (1+sum(sum(repmat((lnm_dist(:,i)),1,Perms)...
                    <=repmat((lnm_dist_null(:,i)),1,Trials)')))/(Trials*Perms+1); 
        elseif strcmp(test,'negative')
            p{n}(i)=...
                    (1+sum(sum(repmat(-1*(lnm_dist(:,i)),1,Perms)...
                    <=repmat(-1*(lnm_dist_null(:,i)),1,Trials)')))/(Trials*Perms+1); 
        end
    end

    %Compute a corrected LNM by z-scoring observed LNM against null dist
    lnm_map{n}=(mean(lnm_dist)-mean(lnm_dist_null))./std(lnm_dist_null);

    %FDR across region p-values
    ind{n}=do_fdr(p{n}',0.05);
    if ~isempty(ind{n})
        fprintf('%d regions significant\n',length(ind{n})); 
    else
        fprintf('Not significant\n');
    end
    fprintf('******\n\n');
end

sig_regions=cell(length(L_set),1); 
sig_regions_name=cell(length(L_set),1);
lesions=cell(length(L_set),1); 
for n=1:length(L_set)
    if ~isempty(ind{n}) %if signifcant regions were found
        %Binary vector with 1 for significant regions
        %Order of regions same as connectivity matrix
        sig_regions{n}=zeros(Nc,1);
        sig_regions{n}(ind{n})=lnm_map{n}(ind{n});

        %List of names of significant regions
        sig_regions_name{n}=Croi(ind{n}); 
        
        %Generate binary image in which voxels comprising significant
        %regions are 1, others 0
        Fname=[L{L_set(n)},'_sig_rois.nii']; 
        vec2nii(sig_regions{n},Fname,atlas_vol,atlas_vals,C2atlas); 

        %Lesion prevalence map 
        x=all_lesions_atlas{L_set(n)}; 
        x=x(sum(x,2)~=0,:); %Remove rows with no lesions
        x_norm=x./repmat(Nc,1,length(C)); 
        lesions{n}=mean(x_norm); 
        lesions{n}=lesions{n}/max(lesions{n}); %normalize to [0,1] for viz
        Fname=[L{L_set(n)},'_lesions.nii']; 
        vec2nii(lesions{n},Fname,atlas_vol,atlas_vals,C2atlas); 

        %Unthresholded and uncorrected LNM map
        %lnm_map{n}=-1*log10(p{n}).*sign(lnm_store(n,:));
        Fname=[L{L_set(n)},'_lnm_map.nii']; 
        vec2nii(lnm_map{n},Fname,atlas_vol,atlas_vals,C2atlas); 
    end
end
