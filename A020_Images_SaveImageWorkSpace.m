%Images_SaveImageWorkSpace

close all;

expno='001_DnaN_TUS_dif_30122014_DnaNsignal_P';

initval=A001_Images_Set_Experiment(expno); %define your paths and files

initval.ImagesWorkspaceName=strcat(initval.basepath,initval.outname,'_Images',num2str(initval.maxfile),'.mat');
[aa,ff,drift]=Get_all_data(initval);
save(initval.ImagesWorkspaceName, 'initval','aa','ff','drift');
disp('done'); 