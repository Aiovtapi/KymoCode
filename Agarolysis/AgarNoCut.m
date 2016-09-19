clc
clear all

%%
init.datapath = 'C:\Users\water\Documents\GitHub\Data\Target Data\Agar 2\';
init.OSslash = '\';
init.flimg = 'stack.tif';

init.stackpath = strcat(init.datapath,init.OSslash,'stack',init.OSslash);
init.meshpath = strcat(init.datapath,'BFC1.mat');

init.pcresize = 1;
init.pctrans = [0,0];

init.Extrabound = 4;
init.strelval = 6;
init.Writebac = 0;
init.TigerCutSR = 1;

%%
init.bacpath = init.datapath;
init.flimgname = {'bacpics.tif'};

javaaddpath('mij.jar')
javaaddpath('ij.jar')

MIJ.start
MIJ.run('Image Sequence...',strcat('open=[',init.stackpath,']'))
MIJ.run('Tiff...',strcat('path=[',init.datapath,init.flimg,']'))
MIJ.closeAllWindows
MIJ.exit

%%
flimg = readtimeseries(strcat(init.datapath,init.flimg));

Ouftiout = load(init.meshpath,'cellList');
Meshdata = Ouftiout.cellList.meshData;
chan = 1;

%%

[Bettermesh,BCellbox,Bacsize,Bacmask,CBacmask,Bacpics,NMBacpics,nflimg] = TigerCut(init,chan,Meshdata,flimg);

%%

frames = size(nflimg,3);
nimgpath = strcat(init.datapath,'TigerCut_',init.flimg);

fprintf('\nWriting frame ')

for frami = 1:frames;
    imwrite(nflimg(:,:,frami),nimgpath, 'writemode', 'append');
    
    if frami>1
        for j=0:log10(frami-1)
            fprintf('\b'); % delete previous counter display
        end
    end
    fprintf('%d', frami);
end