function twopoints=Processing_Pick_Channel(im1,initval,message)

%pick a point in the brightfield%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    imx=max(max(im1)); imn=min(min(im1));
    dipshow(im1,[imn imx]);   hold on 
    disp(message); %show first im
    [r,c]=ginput(1);     %pick a ROI
    dc=initval.kymolength*cos(initval.kymoangle/180*pi)
    dr=initval.kymolength*sin(initval.kymoangle/180*pi)
    r1=r+dr;
    c1=c+dc;
    twopoints=[[r c];[r1 c1]];
    plot([r r1],[c c1],'-o'); %show a line (length, direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%