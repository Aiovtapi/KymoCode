 function [ycor,xcor]=Fine_Tune_ChannelMap(pic,initval);
 %This function finds 'adjustment coordinates' for a a channel picture

 
 %close(gcf);
 figure; pcolor(pic); colormap bone; shading flat;  title('ori picture');
 
 figure;
 edgeroi=pic(:,1:20)-mean(mean(pic(:,1:20))); 
 prf=(mean(edgeroi)).^2; 
 
 %1)set start of map just beyond channel entrance using parabolic fit on edge
 %signal:
 xcor=subpix_step(prf')+initval.entranceoffset;  
 channelroi=pic(:,20:50)';
 
 %2)find cross-channel adjustment using Center-of-mass on x-section
 prf2=mean((channelroi-mean(mean((channelroi))))).^2;
 [~,ycor,~]=Get_1DCOM(prf2); 

 
 dum=1;
 
 