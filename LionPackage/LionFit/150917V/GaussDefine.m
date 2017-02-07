% Define your experiment

exp='RoySim';

if ~exist('Cell')
    Cell=1;
end

switch exp
    
    case 'Old',
        
        %GaussFit
        Bac=num2str(Cell);
        BacStr='Fluo0Chan01Bac0003';
        Mainfolder=strcat(initval.basepath,'Stacks/Tus/');
        Stackpth=strcat(Mainfolder,BacStr);
        d1{1}=readtimeseries(strcat(Stackpth,'/',BacStr,'Im'),'tif'); %read zstack
        
        %GaussCalcs
        MainPathTus=strcat(initval.basepath,'Stacks/Tus/DataMULTI/');
        MainPathdif=strcat(initval.basepath,'Stacks/dif/DataMULTI/');
        
    case 'Mark',
        
        %GaussFit
        folder = {'01','02', '07', '08', '13', '14', '15', '16', '19','20','21','22','27','28','35','36','37','38','39','40','41','42','47','56'};
        Mainfolder='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\Figures\BacPics\';
        Stackpth=strcat('Fluo0Chan01Bac00',folder{Cell},'\');
        Channel=strcat('Fluo0Chan01Bac00',folder{Cell},'Im');
        d1{1}=readtimeseries(strcat(Mainfolder,Stackpth,Channel,'.tif'),'tif');
        
        %GaussCalc
        MainPathTus='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\Figures\BacPics\';
        MainPathdif='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\Figures\BacPics\';
    
    case 'RoySim'
        
        %GaussFit
        Mainfolder='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/';
        Stackpth='1/';
        Channel=num2str(Cell);
        d1{1}=readtimeseries(strcat(Mainfolder,Stackpth,Channel,'.tif'),'tif');
        
        %GaussCalcs
        MainPathTus='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/Results/';
        MainPathdif='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/Results/';
end
