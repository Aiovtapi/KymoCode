function Pixelreshaped = Rolling_ball(array, radius, light_background,smoothing)
    
    % Calculates and subtracts background from array.
    % Arguments:
    % array - uint8 numpy array representing image
    % radius - radius of the rolling ball creating the background
    % light_background - Does image have light background
    % smoothing - Whether the image should be smoothed before creating the background.
    
    if ~exist('light_background')
        light_background = 1;
    end
    if ~exist('smoothing')
        smoothing = 1;
    end
    
    invert = 0;
    if light_background
        invert = 1;
    end


    if radius <= 10
        ball.shrink_factor = 1;
        arc_trim_per = 24;
    elseif radius <= 30
        ball.shrink_factor = 2;
        arc_trim_per = 24;
    elseif radius <= 100
        ball.shrink_factor = 4;
        arc_trim_per = 32;
    else
        ball.shrink_factor = 8;
        arc_trim_per = 40;
    end
    ball = buildball(ball, radius, arc_trim_per);

        
        
    float_array = array;
    float_array = rolling_ball_float_background(float_array,radius,invert,smoothing,ball);
    background_pixels = float_array(:);

    if invert;
        offset = 255.5;
    else
        offset = 0.5;
    end
    
    pixels = uint16(array(:));

    for p = 1:length(pixels);
        value = pixels(p) - (background_pixels(p) + 255) + offset;
        if value < 0
            value = 0;
        end
        if value > 255
            value = 255;
        end

        pixels(p) = uint16(value);
    end
    
    Pixelreshaped = reshape(pixels,size(array));
end


function ball = buildball(ball, ball_radius, arc_trim_per)
    small_ball_radius = ball_radius / ball.shrink_factor;
    if small_ball_radius < 1;
        small_ball_radius = 1;
    end
    
    rsquare = small_ball_radius * small_ball_radius;
    xtrim = floor(arc_trim_per * small_ball_radius) / 100;
    half_width = round(small_ball_radius - xtrim);
    ball.width = (2 * half_width) + 1;
    ball.data = zeros(1,ball.width^2);

    p = 1;
    for y = 1:ball.width;
        for x = 1:ball.width;
            xval = x - half_width;
            yval = y - half_width;
            temp = rsquare - (xval * xval) - (yval * yval);

            if temp > 0
                ball.data(p) = sqrt(temp);
            end
            p = p +1;
        end
    end
end
                
             
%% rolling_ball_float_background

function Pixelreshaped = rolling_ball_float_background(float_array, radius, invert, smoothing, ball)

    % Create background for a float image by rolling a ball over the image

    pixels = float_array(:);
    shrink = ball.shrink_factor > 1;

    if invert
        pixels = -pixels;
    end

    if smoothing
        smoothed_pixels = smooth(reshape(pixels,size(float_array)),3);
        pixels = smoothed_pixels(:);
    end

    pixels = roll_ball(ball,reshape(pixels,size(float_array)));

    if invert
        pixels = -pixels;
    end
    
    Pixelreshaped = reshape(pixels,size(float_array));
end

%% smooth

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

%% roll_ball

function Pixelreshaped = roll_ball(ball, array)

    % Rolls a filtering object over an image in order to find the
    % image's smooth continuous background.  For the purpose of explaining this
    % algorithm, imagine that the 2D grayscale image has a third (height)
    % dimension defined by the intensity value at every point in the image.  The
    % center of the filtering object, a patch from the top of a sphere having
    % radius 'radius', is moved along each scan line of the image so that the
    % patch is tangent to the image at one or more points with every other point
    % on the patch below the corresponding (x,y) point of the image.  Any point
    % either on or below the patch during this process is considered part of the
    % background.
    
    [height, width] = size(array);
    pixels = array(:);
    z_ball = ball.data;
    ball_width = ball.width;
    radius = ball_width/2;
    cache = zeros(width*ball_width);

    for y =-radius:(height + radius);
        next_line_to_write_in_cache = mod((y + radius),ball_width);
        next_line_to_read = y + radius;
        
        if next_line_to_read < height
            src = next_line_to_read * width +1;
            dest = next_line_to_write_in_cache * width+1;
            cache(dest:dest + width) = pixels(src:src + width);
            p = next_line_to_read * width +1;
            for x = 1:width;
                pixels(p) = -inf;
                p = p + 1;
            end
        end
        
        y_0 = y - radius;
        
        if y_0 < 1
            y_0 = 1;
        end
        
        y_ball_0 = y_0 - y + radius;
        y_end = y + radius;
        
        if y_end >= height
            y_end = height - 1;
        end
        
        for x = -radius:(width + radius);
            z = inf;
            x_0 = x - radius;
            if x_0 < 1
                x_0 = 1;
            end
            x_ball_0 = x_0 - x + radius;
            x_end = x + radius;
            
            if x_end >= width
                x_end = width - 1;
            end
            
            y_ball = y_ball_0;
            
            for yp = y_0:(y_end);
                cache_pointer = mod(yp,ball_width)*width + x_0;
                bp = x_ball_0 + y_ball * ball_width;
                
                for xp = x_0:(x_end);
                    z_reduced = cache(cache_pointer) - z_ball(bp);
                    if z > z_reduced
                        z = z_reduced;
                    end
                    cache_pointer = cache_pointer + 1;
                    bp = bp + 1;
                end
                y_ball = y_ball + 1;
            end

            y_ball = y_ball_0;
            
            for yp = y_0:(y_end);
                p = x_0 + yp * width;
                bp = x_ball_0 + y_ball * ball_width;
                for xp = x_0:(x_end);
                    z_min = z + z_ball(bp);
                    if pixels(p) < z_min;
                        pixels(p) = z_min;
                    end
                    p = p + 1;
                    bp = bp + 1;
                end
                y_ball = y_ball + 1;
            end
        end
    end
    
    Pixelreshaped = reshape(pixels,size(array));
end





