clc 
clear all

user = 'MarkPC';
init.viewchan = 'CFP';

init = AgarDefine(user,init);


%%

RoicotrascaFL(init.Agarpath,init.datapath,init.flimgname,strcat(init.datapath,init.beamshape),init.flresize,init.fltrans);
flimg = readtimeseries(strcat(init.datapath,'RIT_',init.flimgname));

Ouftiout = load(strcat(init.datapath,init.meshesfile),'cellList');
Meshdata = Ouftiout.cellList.meshData;

[Bettermesh,Bacmask,Cellbox] = TigerCut(Meshdata,flimg,init,2);

Whichcells = 1:10;
GaussFitSimedit_Agarolysis(init,Bacmask,Whichcells);


%%
Lionresultspath = strcat(init.bacpath,init.flimgname,init.OSslash,'Results',init.OSslash);

for celli = Whichcells;
    
    thismatpath = strcat(Lionresultspath,'Cell_',num2str(celli,'%03.0f'));
    load(thismatpath,'x');
    ld = x;
    spots = size(x,2);
    
    for frami = 1:size(Bettermesh,2);
        
        thismesh = Bettermesh{celli,frami};
        thiscellbox = Cellbox{celli,frami};
        
        thismesh(:,1) = thismesh(:,1) - thiscellbox(1);
        thismesh(:,2) = thismesh(:,2) - thiscellbox(3);
        thismesh(:,3) = thismesh(:,3) - thiscellbox(1);
        thismesh(:,4) = thismesh(:,3) - thiscellbox(3);
        
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
        end
    end
    save(thismatpath,'ld','-append')
end

disp('Projected to Mesh')

