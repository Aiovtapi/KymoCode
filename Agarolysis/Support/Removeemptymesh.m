function Bettermesh = Removeemptymesh(Meshdata)
    
    frames = size(Meshdata,2);
    cells = size(Meshdata{1},2);
    
    ncelli = 1;
    
    % Find cells with empty meshes and skip in Bettermesh
    for celli = 1:cells;
        for frami = 1:frames;
            if numel(Meshdata{frami}{celli}.mesh) > 20;
                Bettermesh{ncelli,frami} = double(Meshdata{frami}{celli}.mesh);
                ncelli = ncelli + 1;
            end
        end
    end
end