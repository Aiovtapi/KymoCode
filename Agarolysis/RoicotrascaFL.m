function RoicotrascaFL(jarpath,Imgspath,Imgname,Beampath,Resizeval,Transval)
%% Presets
if nargin < 5;
    Resizeval = 1;
    Transval = [0,0];
end

%% Loading the images
% Agarolysispth = 'C:\Users\water\Documents\GitHub\KymoCode\Agarolysis\';
% backcorpth = strcat(Agarolysispth, 'Backcor\');
% addpath(backcorpth)
if ~exist('jarpath','var')
    jarpath = 'D:\Users\water_000\Documents\GitHub\KymoCode\Agarolysis\';
end
if ~exist('Imgspath','var')
    Imgspath = 'D:\Users\water_000\Documents\GitHub\Data\141230_dnaN_dif_tus\Fluorescence\CFP\';
end
if ~exist('Imgname','var')
    Imgname = 'CFP.tif';
end
if ~exist('Beampath','var')
    Beampath = 'D:\Users\water_000\Documents\GitHub\Data\141230_dnaN_dif_tus\Fluorescence\BeamShape457.tif';
end

Rballedname = strcat('Rballed_',Imgname);

copyfile(strcat(Imgspath,Imgname),strcat(Imgspath,Rballedname));
disp('Initializing...')
pause(3)

%% Rolling ball

% Add ImageJ and MIJI to path
javaaddpath(strcat(jarpath,'mij.jar'))
javaaddpath(strcat(jarpath,'ij.jar'))

MIJ.start
MIJ.run('Open...',strcat('path=[',Imgspath,Rballedname,']'))
MIJ.run('Subtract Background...', 'rolling=10 stack')
MIJ.run('Save','Tiff...')
MIJ.closeAllWindows
MIJ.exit

disp('Rolling ball finished')

%%

flpth = strcat(Imgspath,Rballedname);
FLinfo = imfinfo(flpth);
num_images = numel(FLinfo);
Outpath = strcat(Imgspath,'RIT_',Imgname);

Beamimg = im2double(imread(Beampath));
maxVAlue= max(max(Beamimg));
NewBeamImg = Beamimg./maxVAlue;


for k = 1:num_images;
    disp(['Processing image ',num2str(k),' out of ',num2str(num_images)])
    
    FLimg = im2double(imread(flpth,k,'Info',FLinfo));
    
    %% Illumunation Correction

    RIimg = double(imdivide(FLimg,NewBeamImg));
    
    %% Transformation
    
    RITimg = imtranslate(imresize(RIimg,Resizeval),Transval,'FIllvalues',0);
    
    %% Write to file
    imwrite(uint16(65535*RITimg),Outpath,'WriteMode','append','Compression','none');

end

disp('Done')
end
