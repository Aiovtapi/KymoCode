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