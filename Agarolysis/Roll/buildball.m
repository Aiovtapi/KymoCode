function ball = buildball(ball, ball_radius, arc_trim_per)
    small_ball_radius = ball_radius / ball.shrink_factor;
    if small_ball_radius < 1;
        small_ball_radius = 1;
    end
    
    rsquare = small_ball_radius * small_ball_radius;
    xtrim = floor(arc_trim_per * small_ball_radius) / 100;
    half_width = round(small_ball_radius - xtrim);
    ball.width = (2 * half_width); %+ 1;
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
                