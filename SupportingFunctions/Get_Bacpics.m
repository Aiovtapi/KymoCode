function [Bacpics, MDchans, FLchans, Cells] = Get_Bacpics(user,exp)

    initval=A001_Images_Set_Experiment(user,exp);

    BacMainfolder = strcat(initval.basepath,initval.FiguresFolder,'BacPics',initval.OSslash);

    FLchans = numel(initval.viewchan);
    MDchans = initval.channelno;
    finddif = strcmp(initval.difchan,initval.viewchan);
    viewFLchans = {initval.viewchan{finddif},initval.viewchan{~finddif}};
    fprintf('\nSaving Bacpics to workspace')

    for mdchan = 1:MDchans
        MDChanMainfolder=strcat(BacMainfolder,'Channel_0',num2str(mdchan),initval.OSslash);
        for flchan = 1:FLchans
            FLChanMainfolder=strcat(MDChanMainfolder,viewFLchans{flchan},initval.OSslash);

            thisdir = dir(FLChanMainfolder);
            thisdir = thisdir(arrayfun(@(x) ~strcmp(x.name(1),'.'),thisdir));
            Ncells = numel(thisdir);
            Cells(mdchan)  = Ncells;

            for celli=1:Ncells

                thisbacfolder=strcat('cell',num2str(celli,'%03.0f'));
                bacseriepth=strcat(FLChanMainfolder,thisbacfolder,initval.OSslash);

                dipbac=readtimeseries(strcat(bacseriepth,'.tif'),'tif');
                uintstack = dip_array(dipbac); %turn into uint16 array

                for frami = 1:size(uintstack,3);
                    Bacpics{mdchan,flchan}{celli,frami} = uintstack(:,:,frami);
                end
            end
        end
    end
end