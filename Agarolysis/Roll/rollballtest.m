[height, width] = size(array);
rarray = array';
pixels = rarray(:); %numpy.float32(array.flatten())
z_ball = ball.data;
ball_width = ball.width;
radius = ball_width/2;
cache = zeros(1,width*ball_width);

for y =-radius :(height + radius - 1);
    next_line_to_write_in_cache(y+9) = mod((y + radius),ball_width);
    next_line_to_read(y+9) = y + radius;

    if next_line_to_read(y+9) < height
        src(y+9) = next_line_to_read(y+9) * width + 1;
        dest(y+9) = next_line_to_write_in_cache(y+9) * width + 1;
        cache(dest(y+9):dest(y+9) + width -1) = pixels(src(y+9):src(y+9) + width - 1);
        p = next_line_to_read(y+9) * width +1;

        pixels(p:p+width-1)=-inf; % -float('inf')
    end
    
	y_0(y+9) = y - radius + 1;
        
    if y_0(y+9) < 1
        y_0(y+9) = 1;
    end

    y_ball_0(y+9) = y_0(y+9) - y + radius;
    
    y_end(y+9) = y + radius + 1;
    if y_end(y+9) >= height
        y_end(y+9) = height;
    end
end