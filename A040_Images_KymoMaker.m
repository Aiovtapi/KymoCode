
%KymoMaker

%This code  allows the user to click a region with a blob and
%tracks a vector from it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%NOTES: 
%'dipimage' should be installed and started ('dipstart') to make this
%code work
%import an 'images' database first! such as 'Images800.mat'
%Jacob Kerssemakers 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below, you can define your own paths and files-------------------------

exp='001_DnaN_TUS_dif_30122014_M';


actions.getdatabase=1;          % default=1 
actions.reloadclicks=0;         %default=0 (new analysis)
actions.firsttime=1;
initval=A001_Images_Set_Experiment(exp);

%Get raw data--------------------------
if actions.getdatabase
ImagesWorkspaceName=strcat(initval.basepath,initval.outname,'_Images',num2str(initval.maxfile),'_',initval.viewchannel,'.mat');
load(ImagesWorkspaceName,'aa', 'ff','drift');
end
%--------------------------------------
 
[~,~,le]=size(aa); 
initval.maxfile=le;
close all;



im1=squeeze(double(dip_array(aa(:,:,1))));
if ~initval.correctdrift
    drift=0*drift;
end
fr_drift=drift(1,:);

%Manual click channels - option to re-load
inname=strcat(initval.basepath,initval.outname_usr);  %separate database storing user actions (limits re-clicking)
if actions.reloadclicks
%option al second-channel loading

if 1
    load(inname,'manypoints'); 
else %explicit other click-file!
    inname2=strcat(initval.basepath,'DnaN_TUS_dif_UserInputs');
    load(inname2,'manypoints'); 
end

else
manypoints=Processing_PickManyChannels(im1,initval,'Pick lower left channel')
if ~actions.firsttime
save(inname, 'manypoints', '-append');
else
save(inname, 'manypoints');
end
end


%% prepare 'presets': corrected start position, reference map or curves etc.

for i=1:initval.channelno
   close all 
   endpoints=[[manypoints(i,1) manypoints(i,2)]; [manypoints(i,3) manypoints(i,4)]];
presets.twopoints=endpoints;

presets.type='BF';
presets.adjustxy=1;
presets.showmap=0;
presets.storeref=1;
presets.useref=0;

[~,presets,~]=Get_Channel_Map(im1, fr_drift, initval,presets);  %process brightfield; use to re-adjust clicked point

presets.adjustxy=0;
presets.storeref=1;

[~,presets,~]=Get_Channel_Map(im1, fr_drift, initval,presets);  % again process brightfield (show&store adjusted result);

%% now choose kymograph channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
presets.type='BF';
% presets.twopoints=Processing_Pick_Channel(im1,initval,'Choose target channel');

presets.adjustxy=1;
presets.storeref=0;
presets.useref=0;
[~, presets,~]=Get_Channel_Map(im1, fr_drift, initval,presets);  %process brightfield; re-adjust clicked point

presets.adjustxy=0;
presets.useref=1;
[~]=Get_Channel_Map(im1, fr_drift, initval,presets);   % again process brightfield (show adjusted result);;

%% Finally, build kymographs using  prior pre-sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
presets.showmap=0;
presets.useref=1;
[kymo_FL,kymo_BF,chanstk_BF,chanstk_FL]= Build_Kymos(aa,ff,drift,initval,presets);

close all;
kymo_BF=fliplr(kymo_BF);
kymo_FL=fliplr(kymo_FL);
figure; pcolor(kymo_BF); shading flat; colormap hot;
figure; pcolor(kymo_FL); shading flat; colormap hot;

%Savefiles----------------------------------------------------------------
initval.WorkspaceOutName=strcat('Exp',exp,'Chan_x',num2str(ceil(presets.twopoints(1,1))),initval.viewchannel,'.mat'); %channel data base
lbl1=strcat(initval.basepath,initval.WorkspaceOutName);  %path+channel database
save(lbl1, 'endpoints', 'presets' ,'initval', 'kymo_FL','kymo_BF','chanstk_BF','chanstk_FL');

%Here we first verify whether we need to create the folder
FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchannel));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchannel))
end

mkdir(strcat(initval.basepath,'/Kymographs/',initval.viewchannel));
lbl3=strcat(initval.basepath,'/Kymographs/',initval.viewchannel,'/Kymograph_FL',initval.WorkspaceOutName(1:end-4),'.tif'); %kymograph
lbl4=strcat(initval.basepath,'/Kymographs/',initval.viewchannel,'/Kymograph_BF',initval.WorkspaceOutName(1:end-4),'.tif'); %kymograph
%For writing Kymos to the Figures folder.
lbl3_Fig=strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchannel,'/Kymograph_FL',initval.WorkspaceOutName(1:end-4),'.tif'); %kymograph
lbl4_Fig=strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchannel,'/Kymograph_BF',initval.WorkspaceOutName(1:end-4),'.tif'); %kymograph

kymim1 = uint8(round(kymo_FL/max(kymo_FL(:))*255 - 1));
imwrite(kymim1,lbl3,'tif');
kymim1 = uint8(round(kymo_BF/max(kymo_BF(:))*255 - 1));
imwrite(kymim1,lbl4,'tif');


imwrite(kymim1,lbl3_Fig,'tif');
imwrite(kymim1,lbl4_Fig,'tif');
end




