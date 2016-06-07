function [BCellbox,bacsize,Bettermesh] = Removeoutbound(BCellbox,bacsize,Bettermesh,flimgsize,frames)
    
    % Find cells whose bounds are outside the FL image
    for frami = 1:frames;
        xlow = find(BCellbox(:,frami,1) < 0);
        xhigh  = find(BCellbox(:,frami,2) > flimgsize(1));
        ylow = find(BCellbox(:,frami,3) < 0);
        yhigh = find(BCellbox(:,frami,4) > flimgsize(2));
    end
    
    % Find cells that are ouf of bound
    fcells = unique([xlow; xhigh; ylow; yhigh])';

    % Remove cells our of bounds of FL image
    for fcelli = fliplr(fcells);
        BCellbox(fcelli,:,:) = [];
        Bettermesh(fcelli,:) = [];
        bacsize(fcelli,:) = [];
    end
end