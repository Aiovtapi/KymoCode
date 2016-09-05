function A020_Images_SaveImageWorkSpace(user,exp)

if nargin<2
    exp='Exp001_DnaN_TUS_dif_01092016_M';
end
if nargin<1
    user = 'MarkPC';
end

initval=A001_Images_Set_Experiment(user,exp); %define your paths and files

ColourNum=size(initval.viewchan,2);

for i=1:ColourNum
initval.ImagesWorkspaceName=strcat(initval.basepath,'Images',num2str(initval.maxfile),'_',initval.outname{i},'.mat');
[aa,ff,drift]=Get_all_data(initval,i);
save(initval.ImagesWorkspaceName, 'initval','aa','ff','drift');

disp(strcat(initval.viewchan{i},' done')); 
end

disp('A20 done')

