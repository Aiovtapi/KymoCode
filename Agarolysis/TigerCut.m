function [Bettermesh,Cellbox] = TigerCut(Meshdata,flimg,init,Extrabound)

disp(sprintf('-----\nOperating TigerCut'))

frames = size(Meshdata,2);
cells = size(Meshdata{1},2);
flimgsize = size(flimg);

[maxX,minX,maxY,minY]=deal(zeros(1,frames));
Bettermesh = cell(cells,frames);
Cellbox = zeros(cells,4);

bacfolder = strcat(init.datapath,'Bacpics',init.OSslash,init.flimgname);

if ~exist(bacfolder,'dir')
    mkdir(bacfolder)
end

for celli = 1:cells;
    for frami = 1:frames;
        Bettermesh{celli,frami} = Meshdata{frami}{celli}.mesh;

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
        maxmesh = max(Bettermesh{celli,frami});
        minmesh = min(Bettermesh{celli,frami});

        maxX(frami) = max(maxmesh(1),maxmesh(3));
        minX(frami) = min(minmesh(1),minmesh(3));
        maxY(frami) = max(maxmesh(2),maxmesh(4));
        minY(frami) = min(minmesh(2),minmesh(4));
    end
    Cellbox(celli,1) = min(minX(:));
    Cellbox(celli,2) = max(maxX(:));
    Cellbox(celli,3) = min(minY(:));
    Cellbox(celli,4) = max(maxY(:));
end


Cellbox(:,1) = Cellbox(:,1) - Extrabound;
Cellbox(:,2) = Cellbox(:,2) + Extrabound;
Cellbox(:,3) = Cellbox(:,3) - Extrabound;
Cellbox(:,4) = Cellbox(:,4) + Extrabound;

RCellbox = round(Cellbox);

% Find cells whose bounds are outside the FL image
fcells = find(RCellbox(:,1) < 0 | RCellbox(:,2) > flimgsize(1) | ...
    RCellbox(:,3) < 0 | RCellbox(:,4) > flimgsize(2));

% Remove cells our of bounds of FL image
for fcelli = fcells;
    Cellbox(fcelli,:) = [];
    RCellbox(fcelli,:) = [];
    Bettermesh(fcelli,:) = [];
end

ncells = size(Bettermesh,1);

disp('Creating Bacpics')

for celli = 1:ncells;
    bacpath=strcat(bacfolder,init.OSslash,'Cell_',num2str(celli,'%03.0f'),'.tif');
    masksize = round([Cellbox(celli,2)-Cellbox(celli,1),Cellbox(celli,4)-Cellbox(celli,3)]/init.pcresize);

    for frami = 1:frames;
        
        % Create mask from mesh
        thismesh = Bettermesh{celli,frames};
        thismesh = [thismesh(:,1:2);thismesh(:,3:4)];
        thiscropmesh = round([thismesh(:,1)-Cellbox(celli,1), thismesh(:,2)-Cellbox(celli,3)]/init.pcresize);
        mask = poly2mask(thiscropmesh(:,2)',thiscropmesh(:,1)',masksize(1),masksize(2));
        
        % Create bacpic from Cellbox
        imageframe = double(flimg(:,:,frami-1));
        croppedimg = imageframe(RCellbox(celli,1):RCellbox(celli,2),RCellbox(celli,3):RCellbox(celli,4));

        % Remove non-cell pixels by applying mask to bacpic
        cimgsize = size(croppedimg);
        nmask = imresize(mask,cimgsize);
        bacpic = uint16(double(nmask).*croppedimg);

        imwrite(bacpic,bacpath,'WriteMode','append','Compression','none')
    end
end

disp(sprintf('TigerCut done \n-----'))


