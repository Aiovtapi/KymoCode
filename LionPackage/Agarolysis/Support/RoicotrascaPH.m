function RoicotrascaPH(jarpath,Imgspath,Imgname,Rval,Tval)
%% Presets
if nargin < 4;
    Rval = 1;
    Tval = [0,0];
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

% Add ImageJ and MIJI to path
javaaddpath(strcat(jarpath,'mij.jar'))
javaaddpath(strcat(jarpath,'ij.jar'))

%% Transformation

imgpath = strcat(Imgspath,Imgname);
imginfo = imfinfo(imgpath);
num_images = numel(imginfo);
Tpath = strcat(Imgspath,'Trans_',Imgname);

RYval = floor(Rval*imginfo(1).Width);
RXval = floor(Rval*imginfo(1).Height);

for k = 1:num_images;
    disp(['Processing image ',num2str(k),' out of ',num2str(num_images)])
    
    img = im2double(imread(imgpath,k,'Info',imginfo));
    
    
    %% Transformation
    
    Timg = img;
    if ~isequal(Rval,1)
        Timg = imresize(Timg,[RXval,RYval],'bilinear');
    end
    
    if ~isequal(Tval,[0,0])
        Timg = imtranslate(Timg,Tval,'FillValues',0);
    end  
    
    % Write to file
    imwrite(uint16(65535*Timg),Tpath,'WriteMode','append','Compression','none');
    
end

% Xstr = num2str(floor(Rval*imginfo(1).Width));
% Ystr = num2str(floor(Rval*imginfo(1).Height));
% Rstr = num2str(Rval);
% 
% MIJ.start
% MIJ.run('Open...',strcat('path=[',imgpath,']'))
% if ~isequal(Rval,0);
%     MIJ.run('Scale...', ['x=',Rstr,' y=',Rstr,' width=',Xstr,' height=',Ystr,' interpolation=',...
%         'Bilinear average create title=',TRname]);
% end
% if ~isequal(Tval,[0,0])
%     MIJ.run('Translate...', ['x=',num2str(Tval(1)),' y=',num2str(Tval(2)),' interpolation=None']);
% end
% MIJ.run('Save','Tiff...')
% MIJ.closeAllWindows
% MIJ.exit

disp('Done')
end
