%Images_SaveImageWorkSpace

close all;

expno='001_DnaN_TUS_dif_30122014_M';

initval=A001_Images_Set_Experiment(expno); %define your paths and files

ColourNum=size(initval.viewchan,2);

for i=1:ColourNum
initval.ImagesWorkspaceName=strcat(initval.basepath,'Images',num2str(initval.maxfile),'_',initval.outname{i},'.mat');
[aa,ff,drift]=Get_all_data(initval,i);
save(initval.ImagesWorkspaceName, 'initval','aa','ff','drift');

disp(strcat(initval.viewchan{i},' done')); 
end

