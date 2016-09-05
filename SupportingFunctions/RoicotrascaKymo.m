 function RoicotrascaKymo(initval,Imgname,Beampath)
%% Presets

OSslash = initval.OSslash;
jarpath = initval.Supportpath;
Imgspath = strcat(initval.basepath,initval.FLpath);
Rval = initval.flresize;
Tval = initval.fltrans;
Raw = initval.rawfolder;

Outpath = strcat(Imgspath,Imgname,OSslash);

if ~exist(Outpath,'dir')
    mkdir(Outpath);
end

% Add ImageJ and MIJI to path
javaaddpath(strcat(jarpath,'mij.jar'))
javaaddpath(strcat(jarpath,'ij.jar'))

%% Rolling ball
Rballedname = strcat('Rballed_',Imgname,'.tif');
% copyfile(strcat(Imgspath,Imgname),strcat(Outpath,Rballedname));
% disp('Copying file...')
pause(3)

MIJ.start
disp('Opening image...')
MIJ.run('Image Sequence...',strcat('open=[',Imgspath,Imgname,OSslash,Raw,']'))
MIJ.run('Tiff...',strcat('path=[',Outpath,Imgname,'.tif]'))

disp('Performing rolling ball background subtraction...')
MIJ.run('Subtract Background...', 'rolling=10 stack')
MIJ.run('Tiff...',strcat('path=[',Outpath,Rballedname,']'))
MIJ.closeAllWindows
MIJ.exit

disp('Rolling ball finished')

%%

imgpath = strcat(Outpath,Rballedname);
imginfo = imfinfo(imgpath);
num_images = numel(imginfo);
RITpath = strcat(Outpath,'RIT_',Imgname,'.tif');

Beamimg = im2double(imread(Beampath));
maxVAlue= max(max(Beamimg));
NewBeamImg = Beamimg./maxVAlue;

RYval = floor(Rval*imginfo(1).Width);
RXval = floor(Rval*imginfo(1).Height);

disp('Performing illumination correction and transformation')
for k = 1:num_images;
    
    img = im2double(imread(imgpath,k,'Info',imginfo));
    
    % Illumunation Correction
    RIimg = double(imdivide(img,NewBeamImg));
    
    %% Transformation
    
    RITimg = RIimg;
    if ~isequal(Rval,1)
        RITimg = imresize(RITimg,[RXval,RYval],'bilinear');
    end
    
    if ~isequal(Tval,[0,0])
        RITimg = imtranslate(RITimg,Tval,'FillValues',0);
    end  
    
    % Write to file
    imwrite(uint16(65535*RITimg),RITpath,'WriteMode','append','Compression','none');
end

disp('Roicotrasca done')
 end

