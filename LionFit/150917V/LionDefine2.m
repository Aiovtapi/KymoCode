function [MainPathTus,MainPathdif] = LionDefine2(exp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Define your experiment



if ~exist('Cell')
    Cell=1;
end

switch exp
    
    case 'Old',
        
        %GaussCalcs
        MainPathTus=strcat(initval.basepath,'Stacks/Tus/DataMULTI/');
        MainPathdif=strcat(initval.basepath,'Stacks/dif/DataMULTI/');
        
    case 'Mark',
        
        %GaussCalc
        MainPathTus='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\Figures\BacPics\';
        MainPathdif='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\Figures\BacPics\';
    
    case 'RoySim'
        
        %GaussFit
        Mainfolder='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/';
        Stackpth='1/';
        Channel=num2str(Cell);
        
        %GaussCalcs
        MainPathTus='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/Results/';
        MainPathdif='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/Results/';
end


end

