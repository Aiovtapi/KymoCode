function [prf,presets,pic]=Get_Channel_Map(im1,fr_drift, initval,presets)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function creates picture based on grid between two points
    %Jacob Kers 2012 
        
    dr=fr_drift(1);     dc=fr_drift(2);                   %drift correction
    r=presets.twopoints(1,1)+dr; c=presets.twopoints(1,2)+dc;
    r1=presets.twopoints(2,1)+dr; c1=presets.twopoints(2,2)+dc;
  
    hw=initval.kymohwidth;
    kl=initval.kymolength;
    radz=initval.kymoangle*pi/180;
    rws= linspace(r, r1,kl);        %create a center line
    cls= linspace(c, c1,kl);
   
    %build grid coordinates,
  % option 2: direct plus adjustments%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   pic=zeros(2*hw+1,kl); 
  rowcoords=zeros(2*hw+1,kl); 
  colcoords=zeros(2*hw+1,kl);
  
  for i=1:kl
        xpos=cls(i);
        ypos=rws(i);
        x_intpol=xpos+cos(radz+pi/2)*((1:2*hw+1)'-hw-1);
        y_intpol=ypos+sin(radz+pi/2)*((1:2*hw+1)'-hw-1);   
        colcoords(:,i)=x_intpol;
        rowcoords(:,i)=y_intpol;
  end
  pic=interp2(im1,rowcoords,colcoords,'linear',0);
%   if presets.type=='FL', 
%       pcolor(pic);  shading flat; colormap hot; 
%   end
  if presets.adjustxy==1
      [rowcor,colcor]=Fine_Tune_ChannelMap(pic,initval);
      presets.twopoints=Adjust_twopoints(rowcor,colcor,presets.twopoints,radz);  
  end
  
  if presets.storeref==1; 
      [r,c]=size(pic);
      pic=repmat(mean(pic')',1,c);  %average along x-direction
      presets.refpic=pic;     
  end;
  
   
  %process picture into profile, depending on channel type
  switch presets.type
  case 'BF'  %BrightField, requires extra enhancement steps-------------
       %1)----- the image is subtracted from the reference picture
        if presets.useref==1;     
          pic=(pic-presets.refpic)./presets.refpic; 
          %pic=(pic)./presets.refpic; 
        end 
        %----------------------------------
        %2) the upper and lower lines of the original imageare taken as background; 
        %they are tiled in an image of the same size as the original. %This
        %image is subtracted from the original.
        [r,c]=size(pic);
        bcklines=[pic(1:2,:) ; pic(r-1:r,:)];  %outer lines selection
        bck=repmat(mean(bcklines),r,1);      pic=pic-bck;
        
        if presets.showmap==1,
            pcolor(pic); colormap bone; shading flat;  
            title('resampled, background corrected picture'); 
        end
        %------------------------------------------
        %3) The inner lines are used to create a derivative edge picture; 
        %The variation in this standard deviation along the Y-direction (perpendicular to the channel)...
        % is taken as a measure for sharp transitions along the channel direction 
        hw=5;
        lo=ceil(r/2)-hw;
        hi=ceil(r/2)+hw;
        bac_pic=pic(lo:hi,:);  %inner lines selection
        for i=2:kl
            edge_pic(:,i)=((bac_pic(:,i)-bac_pic(:,i-1)).^2); 
        end
        finpic=double(dip_array(smooth(edge_pic,1.5)));     
        prf=(std(finpic)-min(std(finpic))); %profile of standard deviations;
        dum=1;
        %-----------------------------------------------
        
      case 'FL'  %Note that fluorescence profiles are NOT Re-normalized
      [rp,cp]=size(pic);
      prf=max(pic)-min(pic);  % fluorescence 
      if presets.showmap==1,  pcolor(pic); colormap bone; shading flat;  title('fluorescence'); end
  end
  
end
