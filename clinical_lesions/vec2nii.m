function vec2nii(img_vec,Name,atlas_vol,atlas_vals,C2atlas)

img=zeros(size(atlas_vol)); 
Nc=length(img_vec); 
for i=1:Nc
    ind_tmp=find(atlas_vals(C2atlas(i))==atlas_vol); 
    img(ind_tmp)=img_vec(i); 
end
mat2nii(img,['./results/',Name],size(img),32,'MNI152_T1_2mm_brain.nii');