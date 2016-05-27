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
TRname = strcat('Trans_',Imgname);

Xstr = num2str(floor(Rval*imginfo(1).Width));
Ystr = num2str(floor(Rval*imginfo(1).Height));
Rstr = num2str(Rval);

MIJ.start
MIJ.run('Open...',strcat('path=[',imgpath,']'))
MIJ.run('Scale...', ['x=',Rstr,' y=',Rstr,' width=',Xstr,' height=',Ystr,' interpolation=',...
    'Bilinear average create title=',TRname]);
MIJ.run('Translate...', ['x=',num2str(Tval(1)),' y=',num2str(Tval(2)),' interpolation=None']);
MIJ.run('Save','Tiff...')
MIJ.closeAllWindows
MIJ.exit

disp('Done')
end
