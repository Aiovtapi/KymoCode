function [Bettermesh,BCellbox,Bacsize,Bacmask,CBacmask,Bacpics,NMBacpics,nflimg] = TigerCut(init,chan,Meshdata,flimg)
    fprintf('-----\nOperating TigerCut')

    flimgsize = size(flimg);
    Bettermesh = Removeemptymesh(Meshdata);
    cells = size(Bettermesh,1);
    meshframes = size(Bettermesh,2);
    
    if numel(flimgsize) == 3;
        if ~meshframes==flimgsize(3)
            disp('Warning: number of frames in FL images does not correspond with number of frames in meshes')
        end
        frames = flimgsize(3);
        if meshframes == 1;
            NBettermesh = cell(cells,frames);
            for celli = 1:cells;
                NBettermesh(celli,:) = repmat(Bettermesh(celli,1),[1,frames]);
            end
            Bettermesh = NBettermesh;
        end
    else
        frames = 1;
    end
    
    nflimg = uint16(zeros(flimgsize(1),flimgsize(2),frames));


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
    [BCellbox,Bacsize, Bettermesh] = Findbound(Bettermesh,Cellbox,cells,frames,init.Extrabound);

    % Remove cells that move out of the immage
    [BCellbox,Bacsize,Bettermesh] = Removeoutbound(BCellbox,Bacsize,Bettermesh,flimgsize,frames,init.Extrabound);

    ncells = size(Bettermesh,1);
    [Bacmask,CBacmask,Bacpics,NMBacpics] = deal(cell(ncells,frames));

    
    if isfield(init,'Writebac')
        Writebac = init.Writebac;
    else
        Writebac = 1;
    end

    if isfield(init,'TigerCutSR')
        TigerCutSR = init.TigerCutSR;
    else
        TigerCutSR = 0;
    end

    if Writebac == 1 ;
        
        fprintf('\nCreating Bacpics')
        fprintf('\nCell: ')
        
        for celli = 1:ncells;
            bacpath=strcat(bacfolder,init.OSslash,'Cell_',num2str(celli,'%03.0f'),init.OSslash);

            if ~exist(bacpath,'dir')
                mkdir(bacpath)
            end
            
            % Display celli number
            if celli>1
                for j=0:log10(celli-1)
                    fprintf('\b');
                end
            end
            fprintf(num2str(celli))

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
    end
    
    if TigerCutSR == 1;
        
        disp(' for SR...')
        fprintf('\nFrame: ')
        
        for frami = 1:frames
            
            if frami>1
                for j=0:log10(frami-1)
                    fprintf('\b');
                end
            end
            fprintf('%d', frami);
            
            % dip_image stack handeling
            if numel(flimgsize) == 2
                imageframe = double(flimg(:,:));
            else
                imageframe = double(flimg(:,:,frami-1));
            end
            
            totalmask = false(flimgsize(1),flimgsize(2));
            
            for celli = 1:ncells
                thismesh = Bettermesh{celli,frami} - init.Extrabound;
                thisrmesh = round([thismesh(:,1:2);thismesh(:,3:4)]);
                mask = poly2mask(thisrmesh(:,1)',thisrmesh(:,2)',flimgsize(2),flimgsize(1));
                dimask = imdilate(mask,strel('disk',init.strelval));
                totalmask = totalmask | dimask;
            end
 
            nflimg(:,:,frami) = uint16(double(totalmask).*imageframe);

        end
    end
    
	fprintf('\nTigerCut done')
end