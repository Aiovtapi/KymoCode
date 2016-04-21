function A000_IlluminationCorrection

%initval = A001_Set_Experiment_PAmCherry(expName);
% clc;
% clear all;

%warning on;

%%--Import and normalize beam profile image
Beamframe = double((imread('BeamShape515c.tif')));

%Beamframe = double((imread(initval.BeamProfile))); 

% fNameToWriteIllumCorrected = '/Volumes/LittleMonster/PhD_Docs/20140103_D003_DnaN_DelTus/Dataset1/Fluorescence_RLBalled_IllumCorrected/Fluorescence_RLBalled_IllumCorrected.tif';
% fNameRBolled = '/Volumes/LittleMonster/PhD_Docs/20140103_D003_DnaN_DelTus/Dataset1/Fluorescence_RollingBalled/Fluorescence_RLBalled.tif';
% 
% 
% fNameToWriteIllumCorrected = 'D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20140127_DnaX_YPet_AND_mCherry_DnaN\Fluorescence\mCherry_Fluorescence_RLBalled_IllumCorrected\mCherry_Fluorescence_RLBalled_IllumCorrected.tif';
% fNameRBolled = 'D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20140127_DnaX_YPet_AND_mCherry_DnaN\Fluorescence\mCherry_DnaN_Rballed.tif';

fNameDirectory='/Users/rleeuw/Work/Data/160205_BN2384_and_Beam_Profiles/2/';
fName='RB/515-100ms-50mWo.tif';
fNameRBolled = strcat(fNameDirectory,fName);
fNameToWriteIllumCorrected = '/Users/rleeuw/Work/Data/160205_BN2384_and_Beam_Profiles/2/IC/515-100ms-33mWo.tif';

if exist(strcat(fNameDirectory,'IC'),'dir')==0
    mkdir(strcat(fNameDirectory,'IC'));
end

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

disp('done')