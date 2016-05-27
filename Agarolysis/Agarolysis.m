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

Agarpath = strcat(init.kymopath,'Agarolysis',init.OSslash);
addpath(genpath(strcat(Agarpath)));
addpath(strcat(init.kymopath,'LionFit',init.OSslash,'150917V'));

init.bfimgname = 'BF.tif';
init.pcimgname = 'PC.tif';
init.flimgname = '457-100ms-10mWo-300G.tif';
init.refimgname = init.pcimgname;
init.maxfile = 421;

init.pcresize = 0.421;
init.pctrans = [-3,62];

% [refimg,flimg]=Get_data(init);

RoicotrascaPH(Agarpath,init.datapath,init.refimgname,0.421,[-3,62])