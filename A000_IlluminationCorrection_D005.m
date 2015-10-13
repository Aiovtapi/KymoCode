function A000_IlluminationCorrection_D005

% clc;
% clear all;

%warning on;

%%-----------Change the folder names to what you have on your computer
%%--Import and normalize beam profile image
Beamframe = double((imread('/Volumes/LittleMonster/PhD_Docs/D005_BeamProfiles/20131128_illumination patterns/BeamProfile_515nm_GaussSmoothed20R.tif')));

%%-----------Change the folder names to what you have on your computer
fNameToWriteIllumCorrected = '/Volumes/LittleMonster/PhD_Docs/D005_Experiments/20131217_DnaN_Dif/TIFFs_DnaN_Dif_004/Fluorescence_Illumcorrected/RBall_Illum_Corrected_Stack.tif';
fNameRBolled = '/Volumes/LittleMonster/PhD_Docs/D005_Experiments/20131217_DnaN_Dif/TIFFs_DnaN_Dif_004/Fluorescence_Rollingballed/FL_RollingBalled.tif';

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