function Division=Division_DoTracking(Division,Divclicks,kymoprops, kymo_BF);  %fill in track posses

r=kymoprops.width;
LL=kymoprops.duration;
 [~,lp]=size(Divclicks);
 
 for i=1:lp   % %for each cluster
     lp-i
     ThisCluster=Divclicks(i);
     startfr=ThisCluster.PosClick.firstfr; %go to start point
     stopfr=ThisCluster.PosClick.lastfr;
     
     lfr=startfr-stopfr;
     
     leftposses=zeros(lfr,1);
     rightposses=zeros(lfr,1);
     frames=zeros(lfr,1);
     
     lastleft=ThisCluster.PosClick.firstleft;  %starting position based on manual clicking
     lastright=ThisCluster.PosClick.firstright; 
     c=0;
     for fr=startfr:stopfr  %for all frames where spot exists:
         c=c+1;
         thisline=kymo_BF(fr,:);                    %next line (unless last of kymograph)       
         frames(c)=fr;
         wd=6; msk=1;
         %LEFT
         lo=ceil(max(lastleft-wd, 1)); %select  lo and hi coords
         hi=ceil(min(lastleft+wd, r)); 
         SpotSoi=thisline(lo:hi);
         SpotSoi=SpotSoi-min(SpotSoi);

         %Center of mass%%%%%%%%%%%%%%%%%%%%%%
         SpotSoi=GaussMask(SpotSoi,msk);
        [~,comc,~]=Get_1DCOM(SpotSoi); 
           
        lastleft=lastleft+comc;       %update last known position
        leftposses(c)=lastleft;
        
        %RIGHT
         lo=ceil(max(lastright-wd, 1)); %select  lo and hi coords
         hi=ceil(min(lastright+wd, r)); 
         SpotSoi=thisline(lo:hi);
         SpotSoi=SpotSoi-min(SpotSoi);

        SpotSoi=GaussMask(SpotSoi,msk);  
        [~,comc,~]=Get_1DCOM(SpotSoi); 
        lastright=lastright+comc;       %update last known position
        rightposses(c)=lastright;   
        
     end
          Division(i).PosKyTrac.frames=frames;
          Division(i).PosKyTrac.left=leftposses;
          Division(i).PosKyTrac.right=rightposses;
          %[i length(frames) length(leftposses) length(rightposses)]
     dum=1;
 end
%%%%%%%%%%%%%%%
