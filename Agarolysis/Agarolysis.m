clc 
clear all

user = 'Mark';

init = AgarDefine(user);


%% Roicotrasca on fluorescence images

% Read FL images 
RoicotrascaFL(init.Agarpath,init.datapath,init.flimgname,init.beampath,init.flresize,init.fltrans);
flimg = readtimeseries(strcat(init.datapath,'RIT_',init.flimgname));

%% Get oufti results and Tigercut

% Load oufti meshdata and compute bacpics with TigerCut
Ouftiout = load(init.meshpath),'cellList');
Meshdata = Ouftiout.cellList.meshData;
[Bettermesh,BCellbox,Bacmask,CBacmask,Bacpics,NMBacpics] = TigerCut(Meshdata,flimg,init,init.Extrabound);

%% Lionfit
[cells,frames]=size(Bacpics);
GaussFitSimedit_Agarolysis(init,Bacpics,Bacmask,cells,frames);

%% ProjectToMesh

Lionresultspath = strcat(init.bacpath,init.flimgname,init.OSslash,'Results',init.OSslash);

figure
showfigure = 0;
AllLnorm = [];

for celli = 1:cells;
    
    thismatpath = strcat(Lionresultspath,'Cell_',num2str(celli,'%03.0f'));
    load(thismatpath,'x');
    
    bx = Removespots(x,CBacmask(celli,:));
    
    ld = bx;
    spots = size(bx,2);
    
    CellLength = zeros(frames,1);
    Lnorm = [];
    
    for frami = 1:frames;
        
        thismesh = Bettermesh{celli,frami};
        thiscellbox = BCellbox(celli,frami,:);
        
        thismesh(:,1) = thismesh(:,1) - thiscellbox(1);
        thismesh(:,2) = thismesh(:,2) - thiscellbox(3);
        thismesh(:,3) = thismesh(:,3) - thiscellbox(1);
        thismesh(:,4) = thismesh(:,4) - thiscellbox(3);
        
        % Mouseclick to view next back, keyboard press for exit preview
        if frami == 1;
            if showfigure == 0;
                thisfigure = imagesc(Bacpics{celli,1});
                hold on
                plot(thismesh(:,1),thismesh(:,2),'w',thismesh(:,3),thismesh(:,4),'w','LineWidth',2)
                for spoti = 1:length(x)
                    plot(x{spoti}(frami,2),x{spoti}(frami,4),'rx','LineWidth',2)
                end
                for spoti = 1:spots;
                    plot(bx{spoti}(frami,2),bx{spoti}(frami,4),'cx','LineWidth',2)
                end
                hold off
                title('Mouseclick to view next bacpic, keypress to exit')

                showfigure = waitforbuttonpress;
            else
                close all
            end
        end
        
        % Find Cell length
        meshlength = length(thismesh);
        [CellLength(frami),~] = projectToMesh(thismesh(meshlength,1),thismesh(meshlength,2),thismesh);
        
        Lnormsp = zeros(spots,1);
        
        % Find LD values by projecting xy values on mesh
        for spoti = 1:spots;
            spotxy = bx{spoti}(frami,:);
            Xval = spotxy(2);
            Yval = spotxy(4);
            
            [Lval,Dval] = projectToMesh(Xval,Yval,thismesh);
            varval = sqrt(spotxy(3)^2+spotxy(5)^2);
            
            ld{spoti}(frami,2) = Lval;
            ld{spoti}(frami,4) = Dval;
            ld{spoti}(frami,3) = varval;
            ld{spoti}(frami,5) = varval;
            
            Lnormsp(spots) = Lval/CellLength(frami);
        end
        Lnorm = [Lnorm; Lnormsp];
    end
    save(thismatpath,'bx','ld','CellLength','Lnorm','-append')
    AllLnorm = [AllLnorm; Lnorm];
end

clear celli CellLength Dval frami ld Lionresultspath Lval Meshlength showfigure spoti spotxy varval x Xval Yval Meshdata thismesh
disp('Projected to Mesh')

