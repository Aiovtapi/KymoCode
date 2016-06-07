clc 
clear all

user = 'Mark';
init.viewchan = 'RFP';

init = AgarDefine(user,init);


%% Roicotrasca on fluorescence images

% Read FL images 
RoicotrascaFL(init.Agarpath,init.datapath,init.flimgname,strcat(init.datapath,init.beamshape),init.flresize,init.fltrans);
flimg = readtimeseries(strcat(init.datapath,'RIT_',init.flimgname));

%% Get oufti results and Tigercut

% Load oufti meshdata and compute bacpics with TigerCut
Ouftiout = load(strcat(init.datapath,init.meshesfile),'cellList');
Meshdata = Ouftiout.cellList.meshData;
[Bettermesh,BCellbox,Bacmask,Bacpics,NMBacpics] = TigerCut(Meshdata,flimg,init,init.Extrabound);

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
    ld = x;
    spots = size(x,2);
    
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
        if showfigure == 0;
            thisfigure = imagesc(Bacpics{celli,1});
            hold on
            plot(thismesh(:,1),thismesh(:,2),'r',thismesh(:,3),thismesh(:,4),'r','LineWidth',2)
            hold off
            title('Mouseclick to view next bacpic, keypress to exit')

            showfigure = waitforbuttonpress;
        else
            close all
        end
        
        % Find Cell length
        meshlength = length(thismesh);
        [CellLength(frami),~] = projectToMesh(thismesh(meshlength,1),thismesh(meshlength,2),thismesh);
        
        Lnormsp = zeros(spots,1);
        
        % Find LD values by projecting xy values on mesh
        for spoti = 1:spots
            spotxy = x{spoti}(frami,:);
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
    save(thismatpath,'ld','CellLength','Lnorm','-append')
    AllLnorm = [AllLnorm; Lnorm];
end

clear celli CellLength Dval frami ld Lionresultspath Lval Meshlength showfigure spoti spotxy varval x Xval Yval Meshdata thismesh
disp('Projected to Mesh')

