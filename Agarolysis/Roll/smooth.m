function new_array = smooth(array, window)

    if ~exist('window')
        window = 3;
    end
    
    [dx, dy] = size(array);
    edge = floor(window/2);
    
    new_array = zeros(dx,dy);
    for i = 1:dx;
        for j = 1:dy;
            window_array = array(max(i-edge,1):min(i+edge,dx),...
                                 max(j-edge,1):min(j+edge,dy));
            new_array(i,j) = mean(window_array(:));
        end
    end
end