function A000_IlluminationCorrection

%initval = A001_Set_Experiment_PAmCherry(expName);
% clc;
% clear all;

%warning on;

%%--Import and normalize beam profile image
Beamframe = double((imread('BeamShape_515nm.tiff')));

%Beamframe = double((imread(initval.BeamProfile))); 

% fNameToWriteIllumCorrected = '/Volumes/LittleMonster/PhD_Docs/20140103_D003_DnaN_DelTus/Dataset1/Fluorescence_RLBalled_IllumCorrected/Fluorescence_RLBalled_IllumCorrected.tif';
% fNameRBolled = '/Volumes/LittleMonster/PhD_Docs/20140103_D003_DnaN_DelTus/Dataset1/Fluorescence_RollingBalled/Fluorescence_RLBalled.tif';
% 
% 
% fNameToWriteIllumCorrected = 'D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20140127_DnaX_YPet_AND_mCherry_DnaN\Fluorescence\mCherry_Fluorescence_RLBalled_IllumCorrected\mCherry_Fluorescence_RLBalled_IllumCorrected.tif';
% fNameRBolled = 'D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20140127_DnaX_YPet_AND_mCherry_DnaN\Fluorescence\mCherry_DnaN_Rballed.tif';

fNameToWriteIllumCorrected = '/Users/rleeuw/Data/20141121/dnaN_dif_Tus_halfmW_5min_80msExp_002/Fluorescence/YFP/YFP_RICorrected.tif';
fNameRBolled = '/Users/rleeuw/Data/20141121/dnaN_dif_Tus_halfmW_5min_80msExp_002/Fluorescence/YFP/YFP.tif';

% FolderExistence = exist(strcat(initval.BasePath,initval.FL_Path_IllumCorrected));
% if FolderExistence == 0
%     mkdir(strcat(initval.BasePath,initval.FL_Path_IllumCorrected));
% end

% figure;
% imagesc(Beamframe)
% axis equal tight
% box on
% colormap(gray)
% axis xy;
maxVAlue= max(max(Beamframe));

NewBeamImg = Beamframe./maxVAlue;

info = imfinfo(fNameRBolled);
num_images = numel(info)


for k = 1:num_images
    fr = imread(fNameRBolled, k, 'Info', info);
       
    fr = double(imdivide(double(fr),NewBeamImg));
    
    AMIN = 0;
    AMAX = 65535;
    
    imwrite(uint16(65535*mat2gray(fr,[AMIN AMAX])),fNameToWriteIllumCorrected,'WriteMode','append','Compression','none');

end