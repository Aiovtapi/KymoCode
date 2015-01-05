function [profile, xax,x_intpol,y_intpol]=Get_Perpendicular_profile(im,xpos,ypos,ori,hw)
%Get a proper cross-section around a max position; in the perpendicular
%direction of ori
    %radians
    
    x_intpol=xpos+cos(ori+pi/2)*((1:2*hw+1)'-hw-1);
    y_intpol=ypos+sin(ori+pi/2)*((1:2*hw+1)'-hw-1);
    
    xr=round(x_intpol);
    yr=round(y_intpol);
    
   %plot(y_intpol,x_intpol, 'r-'); hold on
   
   profile=interp2(im,y_intpol,x_intpol,'linear',0);
    
%    profile=zeros(2*hw,1);
%    for j=1:2*hw
%         profile(j)=im(xr(j),yr(j));
%    end
%     
   
    profile(profile==0)=mean(profile~=0);                        %padding
    q8=ceil(hw/8);
    mn=min([profile(1:q8)' profile(2*hw-q8:2*hw)']);
    profile=profile;                                            %background substraction
    
    
    xax=(1:length(profile))';                   %axis in pixels (y-coordinate only)
end