function [lionval] = LionDefine(exp,initval)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Define your experiment

lionval.channr=7;
lionval.viewchan='CFP';
lionval.viewbac=1;

if ~exist('exp')
    exp = 'Mark';
end

switch exp
    
%     case 'Old',
%         
%         %GaussFit
%         Bac=num2str(Cell);
%         BacStr='Fluo0Chan01Bac0003';
%         Mainfolder=strcat(initval.basepath,'Stacks/Tus/');
%         Stackpth=strcat(Mainfolder,BacStr);
%         d1{1}=readtimeseries(strcat(Stackpth,'/',BacStr,'Im'),'tif'); %read zstack
%         
%         %GaussCalcs
%         MainPathTus=strcat(initval.basepath,'Stacks/Tus/DataMULTI/');
%         MainPathdif=strcat(initval.basepath,'Stacks/dif/DataMULTI/');
        
    case 'Mark',
        
        %GaussFit
        lionval.cropindx=1;
        lionval.bacstring={'','Cell_'};
        
        lionval.difchan='CFP';
        lionval.bacfolder='C:\Users\water\Documents\GitHub\Data\141230_dnaN_dif_tus\Figures\BacPics\';
        lionval.OSslash='\';
        
        % This should be outside of the case as every string introduced are
        % defined in Kymocode
        
        lionval.chanfolder=strcat(lionval.bacfolder,'Channel_0',num2str(lionval.channr),lionval.OSslash);
        lionval.Mainfolder=strcat(lionval.chanfolder,lionval.viewchan,lionval.OSslash);
        lionval.diffolder=strcat(lionval.chanfolder,lionval.difchan,lionval.OSslash);
        % lionval.Stackpth=strcat('Fluo0Chan0',num2str(lionval.channr),'Bac00',num2str(Cell,'%02.0f'));
        
        %GaussCalc
        lionval.MainPathTus='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\Figures\BacPics\';
        lionval.MainPathdif='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\Figures\BacPics\';
    
    case 'RoySim'
        %GaussFit
        lionval.cropindx=0;
        
        lionval.Mainfolder='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/';
        lionval.Stackpth=strcat(num2str(Cell),'/');
        lionval.Channel=num2str(Cell);
        
        %GaussCalcs
        lionval.MainPathTus='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/Results/';
        lionval.MainPathdif='/Users/rleeuw/Work/DataAnalysis/BlurLab/DiffusionTests/Results/';
        
    case 'Roy_MM_Tus_dif'
        %GaussFit     
        lionval.bacstring={'cell'};
        lionval.cropindx=1;
        
        lionval.difchan='CFP';
        lionval.bacfolder='/Users/rleeuw/Work/Data/160220_Working_Stacks_diftus/';
        lionval.OSslash='/';

        lionval.chanfolder=strcat(lionval.bacfolder,'AllCells',lionval.OSslash);
        lionval.Mainfolder=strcat(lionval.chanfolder,lionval.viewchan,lionval.OSslash);
        lionval.diffolder=strcat(lionval.chanfolder,lionval.difchan,lionval.OSslash);
        
        %GaussCalcs
        lionval.MainPathTus='/Users/rleeuw/Work/Data/160220_Working_Stacks_diftus/AllCells/YFP/Results/';
        lionval.MainPathdif='/Users/rleeuw/Work/Data/160220_Working_Stacks_diftus/AllCells/CFP/Results/';
end


end

