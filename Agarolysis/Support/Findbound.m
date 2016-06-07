function [BCellbox,Bettermesh,bacsize] = Findbound(Cellbox,Bettermesh,cells,frames,bound)

    framesize = zeros(cells,frames,2);
    BCellbox = zeros(cells,frames,4);
    bacsize = zeros(cells,2);
    
    for celli = 1:cells;
        % Find the length/height of cell in each frame
        for frami = 1:frames;
            framesize(celli,frami,1) = Cellbox(celli,frami,2)-Cellbox(celli,frami,1);
            framesize(celli,frami,2) = Cellbox(celli,frami,4)-Cellbox(celli,frami,3);
        end
        
        % Find the size of the bacpic by taking max length/height of cell
        bacsize(celli,1) = max(framesize(celli,:,1)) + 2*bound;
        bacsize(celli,2) = max(framesize(celli,:,2)) + 2*bound;
        
        % Create new Cellbox that is bounded
        for frami = 1:frames;
            BCellbox(celli,frami,1) = (Cellbox(celli,frami,2) + Cellbox(celli,frami,1) - bacsize(celli,1))/2 - bound;
            BCellbox(celli,frami,2) = (Cellbox(celli,frami,2) + Cellbox(celli,frami,1) + bacsize(celli,1))/2 + bound;
            BCellbox(celli,frami,3) = (Cellbox(celli,frami,4) + Cellbox(celli,frami,3) - bacsize(celli,2))/2 - bound;
            BCellbox(celli,frami,4) = (Cellbox(celli,frami,4) + Cellbox(celli,frami,3) + bacsize(celli,2))/2 + bound;
            
            Bettermesh{celli,frami} = Bettermesh{celli,frami} - bound;
        end 
    end
end