%function Pixelreshaped = Rolling_ball(array, radius, light_background,smoothing)
    
    % Calculates and subtracts background from array.
    % Arguments:
    % array - uint8 numpy array representing image
    % radius - radius of the rolling ball creating the background
    % light_background - Does image have light background
    % smoothing - Whether the image should be smoothed before creating the background.
    
    %%
    imgpth = 'C:\Users\water\Documents\GitHub\KymoCode\Agarolysis\testim\515-100ms-50mWo-300G.tif';
    array = uint8(im2double(imread(imgpth))*255);
    radius = 10;
    light_background = 0;
    smoothing = 1;
    
    
    %%
    
    
    
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
    
    pixels = uint8(array(:));

    for p = 1:length(pixels);
        value = pixels(p) - (background_pixels(p) + 255) + offset;
        if value < 0
            value = 0;
        end
        if value > 255
            value = 255;
        end

        pixels(p) = uint8(value);
    end
    
    Pixelreshaped = reshape(pixels,size(array));
%end


%% buildball
             
%% rolling_ball_float_background

%% smooth

%% roll_ball
