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
init.pcimgname = 'PC.tif';
init.flimgname = '457-100ms-10mWo-300G.tif';
init.refimgname = init.pcimgname;

flimg = imread(strcat(init.datapath,init.flimgname));
pcimg = imread(strcat(init.datapath,init.pcimgname));

flimgpath = strcat(init.datapath,'SimulTrans_',init.flimgname);
pcimgpath = strcat(init.datapath,'SimulTrans_',init.pcimgname);
%%
rounds = 15;
shiftval = 1;
bound = 20;

xtshift = 0;
ytshift = 0;
xshift = 0;
yshift = 0;

for i = 1:rounds
    xshift = xshift + randi([-shiftval,shiftval]);
    yshift = yshift + randi([-shiftval,shiftval]);
    
    xtshift = xtshift + xshift;
    ytshift = ytshift + yshift;
    
    if abs(xtshift) > bound
        xtshift = xtshift - xshift;
    end
    if abs(ytshift) > bound
        ytshift = ytshift - yshift;
    end
    
    nflimg = imtranslate(flimg,[xtshift, ytshift]);
    npcimg = imtranslate(pcimg,[xtshift, ytshift]);
    
    cnflimg = nflimg(bound+1:end-bound,bound+1:end-bound);
    cnpcimg = npcimg(bound+1:end-bound,bound+1:end-bound);
    
    imwrite(nflimg,flimgpath,'WriteMode','append','Compression','none');
    imwrite(npcimg,pcimgpath,'WriteMode','append','Compression','none');
    
    disp(num2str(i))    
end
    
    