function [outimg] = A000_IlluminationCorrection

%%--Import and normalize beam profile image
Beamframe = double((imread('BeamShape515c.tif')));

fNameDirectory='/Users/rleeuw/Work/Data/160205_BN2384_and_Beam_Profiles/2/';
fName='RB/515-100ms-50mWo.tif';
fNameRBolled = strcat(fNameDirectory,fName);
fNameToWriteIllumCorrected = '/Users/rleeuw/Work/Data/160205_BN2384_and_Beam_Profiles/2/IC/515-100ms-33mWo.tif';



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