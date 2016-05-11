%Rolling Ball Background Subtraction
%Illumination Correction
%Translation
%Scaling
tic
clc
clear all
close all

%% Loading the images
Agarolysispth = 'D:\Users\water_000\Documents\GitHub\KymoCode\Agarolysis\';
Imgpth = 'D:\Users\water_000\Documents\GitHub\Data\141230_dnaN_dif_tus\Fluorescence\CFP\';
imgfolder = strcat(Agarolysispth,'testim\');

flpth = strcat(Imgpth,'CFP.tif');
% flpth = strcat(imgfolder,'457-100ms-10mWo-300G.tif');
flrb = im2double(imread(strcat(imgfolder,'457-100ms-10mWo-300G_Rballed.tif')));
backcorpth = strcat(Agarolysispth, 'Backcor\');
addpath(backcorpth)


FLinfo = imfinfo(flpth);
BFimg = imread(strcat(imgfolder,'BeamShape457.tif'));
PCimg = imread(strcat(imgfolder,'PC.tif'));

% Convert to double for calculation

Beamframe = double(BFimg);
maxVAlue= max(max(Beamframe));
NewBeamImg = Beamframe./maxVAlue;

fNameToWriteIllumCorrected = strcat(Imgpth,'CFPRIauto.tif');
%% Presets
BCpreset = {5,0.01,'sh'};
num_images = numel(FLinfo);

for k = 1:num_images
    disp(['Processing image ',num2str(k),' out of ',num2str(num_images)]);
    
    FLimg = imread(flpth,k,'Info',FLinfo);
    FLdb = im2double(FLimg);
    imgS = size(FLdb);
    

    %% All lines
%     zlength = imgS(1)*imgS(2);
%     noise1 = backcor(linspace(1,zlength,zlength),FLdb,BCpreset{1},BCpreset{2},BCpreset{3});
%     noise2 = backcor(linspace(1,zlength,zlength),FLdb',BCpreset{1},BCpreset{2},BCpreset{3});
%     noisemat1 = reshape(noise1,imgS);
%     noisemat2 = reshape(noise2,fliplr(imgS))';
%     
%     Rballedimg = FLdb-noisemat1;

    %% Line4Line
    noiseall1 = [];
    noiseall2 = [];
    for j = 1: imgS(1);
        noise1 = backcor(linspace(1,imgS(2),imgS(2)),FLdb(j,:),BCpreset{1},BCpreset{2},BCpreset{3});
        noiseall1 = [noiseall1; noise1];
    end
%     for i = 1:imgS(2);
%         noise2 = backcor(linspace(1,imgS(1),imgS(1)),FLdb(:,i),BCpreset{1},BCpreset{2},BCpreset{3});
%         noiseall2 = [noiseall2; noise2];
%     end
    noisemat1 = reshape(noiseall1,imgS);
%     noisemat2 = reshape(noiseall2,fliplr(imgS))';
    Rballedimg = FLdb-noisemat1;%-noisemat2;

    
    Rballedimg(Rballedimg<0)=0;

    %% Illumunation Correction

    RIimg = double(imdivide(Rballedimg,NewBeamImg));
    
    imwrite(uint16(65535*RIimg),fNameToWriteIllumCorrected,'WriteMode','append','Compression','none');

end

% figure
% hold on
% plot(FLdb(:))
% plot(flrb(:))
% plot(RIimg(:))
% plot(noiseall1+noiseall2)
% % plot(noise1)
% legend('Original','ImageJ','New Alllines','Detected noise')
% hold off
% 
% figure; imagesc(FLimg)
% figure; imagesc(RIimg)

disp('Done')
toc
