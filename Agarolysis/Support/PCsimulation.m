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
rval = 1.5;
shearx = -0.01;
sheary = 0.02;
scalex = 1.003;
scaley = 1.001;

nflimg = flimg;
npcimg = pcimg;

fls = size(flimg);
pcs = size(pcimg);

fll = imref2d(fls);
pcl = imref2d(pcs);


for i = 1:rounds
    xshift = -3/fls(2);
    yshift = -17/fls(1);
    
    tformfl = affine2d([scalex, sheary, 0; shearx, scaley, 0; round(xshift*fls(2)), round(yshift*fls(1)), 1]);
	tformpc = affine2d([scalex, sheary, 0; shearx, scaley, 0; round(xshift*pcs(2)), round(yshift*pcs(1)), 1]);
    
    nflimg = imwarp(nflimg,tformfl,'FillValues',0,'OutputView',fll);
    npcimg = imwarp(npcimg,tformpc,'FillValues',0,'OutputView',pcl);
    
    nflimg = imrotate(nflimg,rval,'bilinear');
    npcimg = imrotate(npcimg,rval,'bilinear');
    
%     cnflimg = nflimg(bound+1:end-bound,bound+1:end-bound);
%     cnpcimg = npcimg(bound+1:end-bound,bound+1:end-bound);
    
    imwrite(nflimg,flimgpath,'WriteMode','append','Compression','none');
    imwrite(npcimg,pcimgpath,'WriteMode','append','Compression','none');
    
    
    disp(num2str(i))    
end
    
    