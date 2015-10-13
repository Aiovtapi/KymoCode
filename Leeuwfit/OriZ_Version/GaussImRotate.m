% Rotate stack of images
clear all 
close all
clc

% OriZ: 10r | 8 9r 11 12r 14r 18r 23 24r 27 28r 29 30r 35 36

initval.basepath='/Users/rleeuw/Work/Data/20150322_OriZ_dif/';

Bac='Fluo0Chan02Bac0030';
MainString=strcat(initval.basepath,'DnaN_frominitiationtime/',Bac,'/',Bac,'Im');
D=readtimeseries(MainString,'tif');
data=dip_array(D);

for i=1:size(data,3)
    data(:,:,i)=imrotate(data(:,:,i),180);
    imwrite(data(:,:,i),strcat(MainString,'R',num2str(i),'.tif'),'tif');
end
disp(strcat(num2str(Bac),' done'))
