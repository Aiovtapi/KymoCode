% Rotate stack of images
clear all 
close all
clc

expno='001_DnaN_TUS_dif_30122014_TUSsignal';
initval=A001_Images_Set_Experiment(expno);

Bac='Fluo0Chan02Bac0042';
MainString=strcat(initval.basepath,'StacksLong/dif/',Bac,'/',Bac,'Im');
D=readtimeseries(MainString);
data=dip_array(D);

for i=1:size(data,3)
    data(:,:,i)=imrotate(data(:,:,i),180);
    imwrite(data(:,:,i),strcat(MainString,'R',num2str(i),'.tif'),'tif');
end
