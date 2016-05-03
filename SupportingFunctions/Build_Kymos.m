function [kymo_FL,kymo_BF,chanstk_BF,chanstk_FL]= Build_Kymos(aa,ff,drift,initval, presets)
%Build kymographs from moviestacks. JacobKers 2012

kymo_FL=zeros(initval.maxfile,initval.kymolength);
kymo_BF=zeros(initval.maxfile,initval.kymolength);

chanstk_BF=zeros(2*initval.kymohwidth+1,initval.kymolength,initval.maxfile);  %in pixels
chanstk_FL=zeros(2*initval.kymohwidth+1,initval.kymolength,initval.maxfile);  %in pixels

for j=0:initval.maxfile-2;
initval.maxfile-j;
fr_drift=drift(j+1,:);

%obtain re-sampled values from raw images BRIGHTFIELD---------------------
presets.type='BF';
im1=squeeze(double(dip_array(aa(:,:,j))));
[pr_BF,~,BF]=Get_Channel_Map(im1,fr_drift, initval,presets);  %process brightfield 
kymo_BF(j+1,:)=pr_BF';
%--------------------------------------------------------------------------

%obtain re-sampled values from raw images FLUORESCENCE---------------------
presets.type='FL';
im2=squeeze(double(dip_array(ff(:,:,j))));
[pr_FL,~,FL]=Get_Channel_Map(im2,fr_drift,initval,presets);  %process fluorescence
kymo_FL(j+1,:)=pr_FL';
%----------------------------------------------------------

chanstk_BF(:,:,j+1)=BF;
chanstk_FL(:,:,j+1)=FL;

%plot menu
if mod(j,10)==0
subplot(2,1,1); pcolor(kymo_FL'); shading flat; colormap hot; title('Kymograph Fluorescence');     
subplot(2,1,2); pcolor(kymo_BF'); shading flat; colormap hot; title('Kymograph BrightField Edges');
end
pause(0.02);
end
end
