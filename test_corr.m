%Generate a pair of smoothed images defined on a circular domain. The two 
% images are called Map1 and Map2.  Introduce spatial correlation between 
% the two images by adding one image to the other. More specifcially, 
%Map2 <- (1-alpha)*Map2 + alpha*Map1, where alpha controls the extent of
%spatial correlation. Note that the spatial autocorrelation of Map2 will be
%different than Map1 after this operation and the null hypothesis test should
%account for this difference in spatial autocorrelation.
%A localized difference between Map1 and Map2 is introduced. The difference
%is confined 20% of the maps total area and resides at the center of the
%disk. The difference area is sampled from a third image called Map3.
%Specfically, Map1(x) <- Map3(x)*-1*effect_size and 
%             Map2(x) <- Map3(x)*effect_size
%Where x is the inner area occupying 20% of the total area of the circular
%domain
%We then show that Map1 and Map2 remain highly correlated, despite a
%localized difference being strikingly observable. 
clear all
close all

N=201;   %grid resolution
R=1;     %disk radius
FWHM=8; %in pixels, 

sigma=FWHM/2.3548;

% Create Cartesian grid
x=linspace(-R,R,N);
y=linspace(-R,R,N);
[X,Y]= meshgrid(x, y);

% Convert to polar coordinates
[theta,r]=cart2pol(X,Y);

% Mask outside the disk
mask=(r<=R);

ind=find(mask); 

img1=randn(N,N); 
img1=imgaussfilt(img1,sigma,'Padding','symmetric');
img1(~mask)=0; 

img2=randn(N,N); 
img2=imgaussfilt(img2,sigma,'Padding','symmetric');
img2(~mask)=0; 

%Introduce correlation between the two maps
alpha=0.8; %alpha=1 perfectly correlated
           %alpha=0 independent
img2=(1-alpha)*img2+alpha*img1; %introduce correlation


rval=corr(img1(ind),img2(ind)); 


Perms=1000; 
pval=permutation_null(sigma,N,mask,ind,Perms,alpha,rval,[],[]); 

fprintf('Image correlation: r=%0.2f, p=%0.3f\n',rval,pval); 

%Now randomise an inner disk
prop=0.2; %proportion to randomize
effect_size=0.6; 
R_inner=sqrt(prop*R^2); 
mask_inner=(r<=R_inner);
ind_inner=find(mask_inner);

img_tmp=randn(N,N);
img_tmp=imgaussfilt(img_tmp,sigma,'Padding','symmetric');
img_tmp(~mask)=0;

img1_diff=img1; img1_diff(ind_inner)=effect_size*img_tmp(ind_inner); 
img2_diff=img2; img2_diff(ind_inner)=-effect_size*img_tmp(ind_inner); 

rval_diff=corr(img1_diff(ind),img2_diff(ind)); 

pval_diff=permutation_null(sigma,N,mask,ind,Perms,alpha,rval_diff,ind_inner,effect_size); 

fprintf('Image correlation with difference added: r=%0.2f, p=%0.3f\n',rval_diff,pval_diff); 

Delta=img1_diff-img2_diff; 

hf=figure; hf.Position=[50,50,1200,300]; hf.Color='w';
subplot(1,5,1);
    h=imagesc(img1);
    set(h,'AlphaData',mask); axis off;  axis equal; title('Img1'); 
subplot(1,5,2);
    h=imagesc(img2);
    set(h,'AlphaData',mask); axis off; axis equal; title('Img2'); 
subplot(1,5,3);
    h=imagesc(img1_diff);
    set(h,'AlphaData',mask); axis off;  axis equal; title('Img1 diff'); 
subplot(1,5,4);
    h=imagesc(img2_diff);
    set(h,'AlphaData',mask); axis off; axis equal; title('Img2 diff'); 
subplot(1,5,5);
    h=imagesc(Delta);
    set(h,'AlphaData',mask); axis off; axis equal; title('Delta'); 