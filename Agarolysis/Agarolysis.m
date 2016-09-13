clc 

init = AgarDefine;

if any(strcmp(init.difchan,init.channels(init.viewchannels)))
    chans = [init.viewchannels(strcmp(init.difchan,init.channels(init.viewchannels))),...
        init.viewchannels(~strcmp(init.difchan,init.channels(init.viewchannels)))];
else
    choice = questdlg('Dif-channel not selected. Cell orientation may not be defined. Do you want to continue?',...
        'Warning','Yes','No','No');
    switch choice
        case 'No'
            return
    end
end

[Bacpics, NMBacpics, Bacmesh] = deal(cell(numel(chans),1)); 


%% Roicotrasca on fluorescence images
for chan = chans
   

    % Read FL images 
    RoicotrascaFL(init.OSslash,init.Agarpath,init.datapath,init.flimgname{chan},...
        init.beampath{chan},init.flresize,init.fltrans);
    flimg = readtimeseries(strcat(init.datapath,'Edited Images',init.OSslash,'RIT_',init.flimgname{chan}));

    %% Get oufti results and Tigercut

    % Load oufti meshdata and compute Bacpics{chan} with TigerCut
    Ouftiout = load(init.meshpath,'cellList');
    Meshdata = Ouftiout.cellList.meshData;
    
    % Use TigerCutSimple if meshes and cellboxes are already calculated
    if chan == 1;
        [Bettermesh,BCellbox,Bacsize,Bacmask,CBacmask,Bacpics{chan},NMBacpics{chan}]...
            = TigerCut(init,chan,Meshdata,flimg);
    elseif chan > 1;
        [Bacpics{chan},NMBacpics{chan}] = TigerCutSimple(init,chan,OBettermesh,OBCellbox,OBacsize,flimg);
    end
    
    %S Save mesh and cellbox for TigerCutSimple
    if chan == chans(1)
        OBettermesh = Bettermesh;
        OBCellbox = BCellbox;
        OBacsize = Bacsize;
    end
    
    %% Lionfit
    [cells,frames]=size(Bacpics{chan});
    
    IPTPvalue(chan) = LionfitcontrolUI(init,Bacpics{chan},Bacmask,cells,chan);
    close all
    
    % Check if rotation information already exists. If not, rotation
    % information will be generated and saved
    if strcmp(init.channels{chan},init.difchan)
        disp('Calculated Bacpic rotation')
        imgflip = GaussFitSimedit_AgarRotate(init,chan,Bacpics{chan},cells,1,IPTPvalue(chan));
        flippedchan = [chans ;zeros(size(chans))];
        save(strcat(init.datapath,'imgflip.mat'),'imgflip','flippedchan')
    else
        if exist(strcat(init.datapath,'imgflip.mat'),'file')
            load(strcat(init.datapath,'imgflip.mat'))
        else
        	choice = questdlg('No imgflip file found, cell orientation not defined. Do you want to continue?',...
                'Warning','Yes','No','No');
            switch choice
                case 'No'
                    return
            end
        end
    end
    
    
    % Rotation of the cells according to the dif-site information
    
    if exist(strcat(init.datapath,'imgflip.mat'),'file')
        % Only flip Bettermesh and BCellbox once
        if chan == chans(1)
            disp('Flipping Meshes...')
            [Bettermesh, BCellbox, Bacsize] = AgarRotate1(init, Bettermesh, cells, frames, imgflip);
        end
        
        % Flipping the bacpics and masks
        disp('Flipping Bacpics...')
        [Bacmask, CBacmask, Bacpics{chan}, NMBacpics{chan}] = ...
            AgarRotate2(init, Bacmask, CBacmask, Bacpics{chan}, NMBacpics{chan},...
            cells, frames, imgflip, chan);
        
        % Save flipped information
        flippedchan(2,chan) = 1;
        save(strcat(init.datapath,'imgflip.mat'),'flippedchan','-append')
    end
       
    DataStruct = GaussFitSimedit_Agarolysis(init,chan,Bacpics{chan},Bacmask,cells,frames,IPTPvalue(chan));

    %% Get indivudual meshes

    for celli = 1:cells
        for frami = 1:frames;
            thismesh = Bettermesh{celli,frami} + init.Extrabound;
            thiscellbox = BCellbox(celli,frami,:);

            Bacmesh{chan}{celli,frami}(:,1) = thismesh(:,1) - thiscellbox(1);
            Bacmesh{chan}{celli,frami}(:,2) = thismesh(:,2) - thiscellbox(3);
            Bacmesh{chan}{celli,frami}(:,3) = thismesh(:,3) - thiscellbox(1);
            Bacmesh{chan}{celli,frami}(:,4) = thismesh(:,4) - thiscellbox(3);
        end
    end
    
    %% ProjectToMesh

    Lionresultspath = strcat(init.bacpath,init.flimgname{chan},init.OSslash,'Results',init.OSslash);
    AllLnorm = [];

    for celli = 1:cells;

%         thismatpath = strcat(Lionresultspath,'Cell_',num2str(celli,'%03.0f'));
%         load(thismatpath,'x');
        
        x = DataStruct(chan,celli).x;
        bx = Removespots(x,CBacmask(celli,:));
        X{chan,celli} = x;
        BX{chan,celli} = bx;
        
        ld = bx;
        spots = size(bx,2);

        CellLength = zeros(frames,1);
        Lnorm = [];
        
        for frami = 1:frames;

            % Find Cell length
            thisbacmesh = Bacmesh{chan}{celli,frami};
            meshlength = length(thisbacmesh);
            [CellLength(frami),~] = projectToMesh(thisbacmesh(meshlength,1),thisbacmesh(meshlength,2),thisbacmesh);

            Lnormsp = [];

            % Find LD values by projecting xy values on mesh
            for spoti = 1:spots;
                spotxy = bx{spoti}(frami,:);
                Xval = spotxy(2);
                Yval = spotxy(4);

                [Lval,Dval] = projectToMesh(Xval,Yval,thisbacmesh);
                varval = sqrt(spotxy(3)^2+spotxy(5)^2);

                ld{spoti}(frami,2) = Lval;
                ld{spoti}(frami,4) = Dval;
                ld{spoti}(frami,3) = varval;
                ld{spoti}(frami,5) = varval;

                Lnormsp = [Lnormsp; Lval/CellLength(frami)];
            end
            Lnorm = [Lnorm; Lnormsp];
        end
        
        DataStruct(chan,celli).bx = bx;
        DataStruct(chan,celli).ld = ld;
        DataStruct(chan,celli).CellLength = CellLength;
        DataStruct(chan,celli).Lnorm = Lnorm;
        
        AllLnorm = [AllLnorm; Lnorm];
        clear x bx
    end

    Lvalnorm{chan} = AllLnorm;
    
    clear celli CellLength Dval frami Lionresultspath Lval Meshlength showfigure spoti...
        spotxy varval Xval Yval Meshdata thismesh ld Lnorm Lnormsp meshlength...
        spots thiscellbox thisfigure thismatpath AllLnorm thisbacmesh
    
    disp('Projected to Mesh')
end

save(strcat(init.datapath,init.OSslash,'Results.mat'),'DataStruct')


%% View bacpics with meshes and spots, find faulty cells

skip = 0;
fcelli = [];
celli = 1;

f = figure('Name','Agarolysis','NumberTitle','off',...
        'Visible','on','Position',[350 350 1200, 400]);
while celli <= cells;
    
    if skip == 0;
        [skip,fault,previous,Rspot] = ViewbacUI2(f,Bacpics,Bacmesh,X,BX,celli,init.flimgname);
    end
    
    % Remove clicked spots, new bx and ld are saved as rbx and rld
	if ~numel(Rspot) == 0
        Removespotsui(init,celli,Rspot)
    end
    
    % Save faulty cells
    if fault == 1;
        fcelli = [fcelli, celli];
    end
    
    % Go to previous cell
    switch previous
        case 0; celli = celli + 1;
        case 1
            if celli > 1
                celli = celli - 1;
            end
    end
    
    if celli == cells + 1;
        done = questdlg('Advance?','Menu','Yes','No, Go back','Yes');
        if strcmp(done,'No, Go back')
            switch skip
                case 0; celli = celli - 1;
                case 1; skip = 0; celli = 1;
            end
        end
    end
end
close(f)

faultycells = unique(fcelli);
fpath = strcat(init.bacpath,'fcells.mat');
fremoved = 0;
save(fpath,'faultycells','fremoved')

%% Remove faulty cells

if ~isempty(faultycells)
    RemoveCells(init,cells,faultycells,fpath);
end
disp('Operation done')

clear chan chans celli done previous skip frames fremoved fault fclii




