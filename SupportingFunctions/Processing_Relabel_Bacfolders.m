function Processing_Relabel_Bacfolders(exp,user)

if nargin<1, exp='001_DnaN_TUS_dif_30122014_M';
end


initval=A001_Images_Set_Experiment(user,exp);

Nnumbers = 3; % Amount of numbers at the end of the folders

bacpath=strcat(initval.basepath,initval.FiguresFolder,'BacPics',initval.OSslash);

for m = 1:initval.channelno;
    
    channelpath=strcat(bacpath,'Channel_',num2str(m,'%02.0f'),initval.OSslash);
    
    for i = 1:size(initval.viewchan,2)
        
        colorpath = strcat(channelpath,initval.viewchan{i},initval.OSslash);
        bacfolders=dir(colorpath);
        Nbacs = size(bacfolders,1);

        foldername = bacfolders(end).name;
        foldername = foldername(1:length(foldername)-Nnumbers);

        k = 1;
        for j = 1:Nbacs;
            if length(bacfolders(j).name) > length(foldername);

                newname = strcat(foldername,num2str(k,strcat('%0',num2str(Nnumbers),'.0f')));

                if ~strcmp(bacfolders(j).name,newname);
                    movefile(strcat(colorpath,bacfolders(j).name),strcat(colorpath,newname));
                end
                k = k+1;
            end
        end 
    end
end
