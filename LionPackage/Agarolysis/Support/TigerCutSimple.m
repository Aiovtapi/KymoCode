function [Bacpics,NMBacpics] = TigerCutSimple(init,chan,Bettermesh,BCellbox,Bacsize,flimg)

    ncells = size(Bettermesh,1);
    frames = size(Bettermesh,2);
    flimgsize = size(flimg);
    
    bacfolder = strcat(init.bacpath,init.flimgname{chan});
    if ~exist(bacfolder,'dir')
        mkdir(bacfolder)
    end
    
    [Bacmask,CBacmask,Bacpics,NMBacpics] = deal(cell(ncells,frames));

    fprintf('\nCreating Bacpics')
    fprintf('\nCell: ')

    for celli = 1:ncells;
        
        % Display celli number
        if celli>1
            for j=0:log10(celli-1)
                fprintf('\b');
            end
        end
        fprintf(num2str(celli))
        
        bacpath=strcat(bacfolder,init.OSslash,'Cell_',num2str(celli,'%03.0f'),init.OSslash);
        
        if ~exist(bacpath,'dir')
            mkdir(bacpath)
        end
        
        for frami = 1:frames;
            
            % dip_image stack handeling
            if numel(flimgsize) == 2
                imageframe = double(flimg(:,:)); % ,frami-1));
            else
                imageframe = double(flimg(:,:,frami-1));
            end
            
            % Set values for current cell and frame
            thismesh = Bettermesh{celli,frami};
            thisBbox = squeeze(BCellbox(celli,frami,:));
            thisbacsize = Bacsize(celli,:);
            
            % create bacpic and save mask
            [mmask, nmask, bacpic,croppedimg] = Createbac(init,imageframe,thismesh,thisBbox,thisbacsize,bacpath,frami);          
            Bacmask{celli,frami} = nmask;
            CBacmask{celli,frami} = mmask;
            Bacpics{celli,frami} = bacpic;
            NMBacpics{celli,frami} = croppedimg;
        end    
    end
    
    fprintf('\nTigerCut done')
end