function [Bettermesh,BCellbox,Bacsize,Bacmask,CBacmask,Bacpics,NMBacpics] = TigerCut(init,chan,Meshdata,flimg)
    disp(sprintf('-----\nOperating TigerCut'))

    Bettermesh = Removeemptymesh(Meshdata);
    frames = size(Bettermesh,2);
    cells = size(Bettermesh,1);
    flimgsize = size(flimg);

    Cellbox = zeros(cells,frames,4);

    bacfolder = strcat(init.bacpath,init.flimgname{chan});

    if ~exist(bacfolder,'dir')
        mkdir(bacfolder)
    end

    for celli = 1:cells;
        for frami = 1:frames;

            % Add translation to meshes
            if ~isequal(init.pcresize,1)
                Bettermesh{celli,frami} = double(init.pcresize*Bettermesh{celli,frami});
            end

            if ~isequal(init.pctrans,[0,0])
                Bettermesh{celli,frami}(:,1) = init.pctrans(1) + Bettermesh{celli,frami}(:,1);
                Bettermesh{celli,frami}(:,2) = -init.pctrans(2) + Bettermesh{celli,frami}(:,2);
                Bettermesh{celli,frami}(:,3) = init.pctrans(1) + Bettermesh{celli,frami}(:,3);
                Bettermesh{celli,frami}(:,4) = -init.pctrans(2) + Bettermesh{celli,frami}(:,4);
            end  

            % Find mesh maxima and minima
            maxmesh = round(max(Bettermesh{celli,frami}));
            minmesh = round(min(Bettermesh{celli,frami}));

            Cellbox(celli,frami,1) = min(minmesh(1),minmesh(3));            
            Cellbox(celli,frami,2) = max(maxmesh(1),maxmesh(3));
            Cellbox(celli,frami,3) = min(minmesh(2),minmesh(4));
            Cellbox(celli,frami,4) = max(maxmesh(2),maxmesh(4));
        end
    end

    % Find size of bacpic and the boundary indeces for each frame
    [BCellbox,Bettermesh,Bacsize] = Findbound(Cellbox,Bettermesh,cells,frames,init.Extrabound);

    % Remove cells that move out of the immage
    [BCellbox,Bacsize,Bettermesh] = Removeoutbound(BCellbox,Bacsize,Bettermesh,flimgsize,frames);

    ncells = size(Bettermesh,1);
    [Bacmask,CBacmask,Bacpics,NMBacpics] = deal(cell(ncells,frames));

    disp('Creating Bacpics')

    for celli = 1:ncells;
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
    
    disp(sprintf('TigerCut done \n-----'))
end