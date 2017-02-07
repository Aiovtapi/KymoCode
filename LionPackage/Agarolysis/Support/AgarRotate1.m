function [Bettermesh, BCellbox, Bacsize] = ...
    AgarRotate1(init, Bettermesh, cells, frames, imgflip)

pcsize = size(imread(strcat(init.datapath,init.PCimgname)))*init.pcresize;
bpcsize = pcsize + 2*init.Extrabound;
Cellbox = zeros(cells,frames,4);

for celli = 1:cells;
    for frami = 1:frames;
        if imgflip(celli) == 1;
        % Find opposite values of Bettermesh and BCellbox

            Bettermesh{celli,frami}(:,1) = pcsize(2) - Bettermesh{celli,frami}(:,1);
            Bettermesh{celli,frami}(:,2) = pcsize(1) - Bettermesh{celli,frami}(:,2);
            Bettermesh{celli,frami}(:,3) = pcsize(2) - Bettermesh{celli,frami}(:,3);
            Bettermesh{celli,frami}(:,4) = pcsize(1) - Bettermesh{celli,frami}(:,4);
        end

%             BCellbox(celli,frami,1) = round(bpcsize(2) - BCellbox(celli,frami,1));
%             BCellbox(celli,frami,2) = round(bpcsize(2) - BCellbox(celli,frami,2));
%             BCellbox(celli,frami,3) = round(bpcsize(1) - BCellbox(celli,frami,3));
%             BCellbox(celli,frami,4) = round(bpcsize(1) - BCellbox(celli,frami,4));

        % Find mesh maxima and minima
        maxmesh = round(max(Bettermesh{celli,frami}));
        minmesh = round(min(Bettermesh{celli,frami}));

        Cellbox(celli,frami,1) = min(minmesh(1),minmesh(3));            
        Cellbox(celli,frami,2) = max(maxmesh(1),maxmesh(3));
        Cellbox(celli,frami,3) = min(minmesh(2),minmesh(4));
        Cellbox(celli,frami,4) = max(maxmesh(2),maxmesh(4));

    end
end
    [BCellbox,Bacsize] = Findbound(Cellbox,cells,frames,init.Extrabound);
end