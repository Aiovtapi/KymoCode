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
    pixels = array(:); %numpy.float32(array.flatten())
    z_ball = ball.data;
    ball_width = ball.width;
    radius = ball_width/2;
    cache = zeros(1,width*ball_width);
    
    % y = x and x = y, width = height and height = width
    for y =-radius :(height + radius - 1);
        next_line_to_write_in_cache = mod((y + radius),ball_width);
        next_line_to_read = y + radius;

        if next_line_to_read < height
            src = next_line_to_read * width + 1;
            dest = next_line_to_write_in_cache * width + 1;
            cache(dest:dest + width -1) = pixels(src:src + width - 1);
            p = next_line_to_read * width +1;
            pixels(p:p+width-1)=-inf; % -float('inf')
        end

        y_0  = y - radius + 1;

        if y_0  < 1
            y_0  = 1;
        end

        y_ball_0  = y_0  - y + radius;

        y_end  = y + radius + 1;
        if y_end >= height
            y_end = height;
        end
        
        for x = -radius:(width + radius -1);
            z = inf; %float('inf')
            x_0 = x - radius +1;
            if x_0 < 1
                x_0 = 1;
            end
            x_ball_0 = x_0 - x + radius;
            
            x_end = x + radius + 1;
            if x_end >= width
                x_end = width;
            end
            
            y_ball = y_ball_0;
            
            %%%%%%%%%%%%
            x
            
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