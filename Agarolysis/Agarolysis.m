clc 
clear all

user = 'MarkPC';

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
init.pcimgname = 'PC.tif';
init.flimgname = '457-100ms-10mWo-300G.tif';
init.refimgname = init.pcimgname;
init.maxfile = 421;

init.pcresize = 0.421;
init.pctrans = [-3,62];

% [refimg,flimg]=Get_data(init);

%%
RoicotrascaPH(init.Agarpath,init.datapath,init.refimgname,0.421,[-3,62]);
RoicotrascaFL(init.Agarpath,init.datapath,'CFP.tif',strcat(init.datapath,'BeamShape457.tif'),0,[0,0]);

%%

flimg = imread(strcat(init.datapath,init.flimgname));
pcimg = readtimeseries(strcat(init.datapath,'Trans_',init.pcimgname));

flimg2 = readtimeseries(strcat(init.datapath,'RIT_CFP.tif'));

if isa(pcimg,'uint8')
    pcimg = uint16(65535*im2double(pcimg));
end

%%

FLimgsize = size(flimg);




