clc 
clear all

user = 'Mark';

switch user
    case 'Mark'
        init.OSslash = '\';
        init.kymopath = 'C:\Users\water\Documents\GitHub\KymoCode\';
        init.datapath = 'C:\Users\water\Documents\GitHub\Data\DnaN_dif_Tus_AgarPad\';
    case 'MarkPC'
        init.OSslash = '\';
        init.kymopath = 'D:\Users\water_000\Documents\GitHub\KymoCode\';
        init.datapath = 'D:\Users\water_000\Documents\GitHub\Data\DnaN_dif_Tus_AgarPad\';
end

init.Agarpath = strcat(init.kymopath,'Agarolysis',init.OSslash);
addpath(genpath(strcat(init.Agarpath)));
addpath(strcat(init.kymopath,'LionFit',init.OSslash,'150917V'));

init.bfimgname = 'BF.tif';
init.pcimgname = 'Simul_PC.tif';
init.flimgname = 'Simul_457-100ms-10mWo-300G.tif';
init.beamshape = 'BeamShape457.tif';
init.meshesfile = 'Simul_PC.mat';
init.refimgname = init.pcimgname;
init.maxfile = 421;

init.pcresize = 0.421;
init.pctrans = [0,0];
init.flresize = 1;
init.fltrans = [2,-63];


%%

RoicotrascaFL(init.Agarpath,init.datapath,init.flimgname,strcat(init.datapath,init.beamshape),init.flresize,init.fltrans);
flimg = readtimeseries(strcat(init.datapath,'RIT_',init.flimgname));

Ouftiout = load(strcat(init.datapath,init.meshesfile),'cellList');
Meshdata = Ouftiout.cellList.meshData;

[Bettermesh,Cellbox] = TigerCut(Meshdata,flimg,init,2);



% RoicotrascaPH(init.Agarpath,init.datapath,init.refimgname,init.pcresize,init.pctrans);
% pcimg = imread(strcat(init.datapath,'Trans_',init.pcimgname));

% if isa(pcimg,'uint8')
%     pcimg = uint16(65535*im2double(pcimg));
% end

% flimg = double(flimg);
% pcimg = double(pcimg);
% 
% minX = min(size(pcimg,1),size(flimg,1));
% minY = min(size(pcimg,2),size(flimg,2));
% 
% if size(pcimg,1)>minX
%     pcimg(minX+1:end,:)=[];
% end
% if size(flimg,1)>minX
%     flimg(minX+1:end,:)=[];
% end
% if size(pcimg,2)>minY
%     pcimg(:,minY+1:end)=[];
% end
% if size(flimg,2)>minY
%     flimg(:,minY+1:end)=[];
% end
% 
% imwrite(uint16(flimg),strcat(init.datapath,'Trans_',init.flimgname));
% imwrite(uint16(pcimg),strcat(init.datapath,'Trans_',init.pcimgname));

%%

% load(strcat(init.datapath,'Trans_PC.mat'),'cellList')
