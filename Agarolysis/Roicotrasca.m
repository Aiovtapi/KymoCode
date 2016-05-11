%Rolling Ball Background Subtraction
%Illumination Correction
%Translation
%Scaling

clc
clear all
close all

%% Loading the images
Agarolysispth = 'C:\Users\water\Documents\GitHub\KymoCode\Agarolysis\';
imgfolder = strcat(Agarolysispth,'testim\');
flpth = strcat(imgfolder,'515-100ms-50mWo-300G.tif');
backcorpth = strcat(Agarolysispth, 'Backcor\');
addpath(backcorpth)


FLinfo = imfinfo(flpth);
BFimg = imread(strcat(imgfolder,'BeamShape515.tif'));
PCimg = imread(strcat(imgfolder,'PC.tif'));

% Convert to double for calculation

Beamframe = double(BFimg);
maxVAlue= max(max(Beamframe));
NewBeamImg = Beamframe./maxVAlue;


%% Presets
BCpreset = {4,0.1,'sh'};
num_images = numel(FLinfo);

for k = 1:num_images
    FLimg = imread(flpth,k,'Info',FLinfo);
    FLdb = im2double(FLimg);
    imgS = size(FLdb);
    

    % Rolling Ball
    zlength = imgS(1)*imgS(2);
    noise = backcor(linspace(1,zlength,zlength),FLdb(:),BCpreset{1},BCpreset{2},BCpreset{3});
    noisemat = reshape(noise,imgS);

    Rballedimg = FLdb-noisemat;
    Rballedimg(Rballedimg<0)=0;

    %% Illumunation Correction

    RIimg = double(imdivide(Rballedimg,NewBeamImg));
    
%     AMIN = 0;
%     AMAX = 65535;
%     imwrite(uint16(65535*mat2gray(fr,[AMIN AMAX])),fNameToWriteIllumCorrected,'WriteMode','append','Compression','none');

end

