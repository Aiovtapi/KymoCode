function manypoints=Processing_PickManyChannels(im1,initval,message)

%pick a row of points point in the brightfield%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
confirm = 0;    
manypoints=[];
imx=max(max(im1)); imn=min(min(im1));

while ~(confirm == 1)

    dipshow(im1,[imn imx]);   hold on 
    disp(message); %show first im
    disp('Right click to restart')
    [r0,c0]=ginput(1);     %pick a ROI
    
    no=initval.channelno;
    dx=initval.channeldistance;
    shifts=dx*[0:1:no-1];    
    for i=1:no
    
        dc=initval.kymolength*cos(initval.kymoangle/180*pi); %longi shift
        dr=initval.kymolength*sin(initval.kymoangle/180*pi);

        dr1=shifts(i)*sin((initval.kymoangle-90+initval.perpadjust)/180*pi);  %perp shift
        dc1=shifts(i)*cos((initval.kymoangle-90+initval.perpadjust)/180*pi);
        
        r1=r0+dr1;
        c1=c0+dc1;

        r2=r0+dr+dr1;
        c2=c0+dc+dc1;
        manypoints(i,:)=[[r1 c1] [r2 c2]];
        plot([r1 r2],[c1 c2],'-o'); %show a line (length, direction)
       
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 dum=1;
 [~,~,confirm]=ginput(1);
 close all
end