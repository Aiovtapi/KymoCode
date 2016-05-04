    function ReplicationCluster=RepliCluster_DoTrackingQuick(ReplicationCluster,kymoprops,kymo_FL,initval);
    %ReplicationCluster
%----------------------------------------------------------------
% name (genealogic number; 2^n indicates mother cell) 
% fate  'disassembled' ; 'exit' ;'nonexistent' etc.
% linkedbac: label of associated division cycle
% PosClick
%     firstframe: frameno
%     lastframe: frameno
%     firstpos: pixel pos
%     lastpos: pixel pos
% PosKyTracCom
%     frames: frameno
%     clickpos: interpolated between start/stop clicks, pixels
%     trackpos: Centroid tracked, pixels
% PosKyTracGauss
%     spot1pos: Gauss 1 position, pixels
%     spot1amp: Gauss 1 peak height, a.u.
%     spot1sig: common Gauss sigma (user set), pixels
%     spot2pos: Gauss 2 position, pixels
%     spot2amp: Gauss 2 peak height, a.u.
%     spot2sig: common Gauss sigma (user set), pixels 
% Pos2DPreTrac
%     X0: Spot 1 X-position, pixels
%     X1: Spot 2 X-position, pixels
%     Y0: Spot 1 Y-position, pixels
%     Y1: Spot 2 Y-position, pixels
%     Bck: Background level
%     Pk1: Spot 1 Peak, a.u.
%     Pk2: Spot 2 Peak, a.u.
%     spots: spot number present, according to rejection criteria
%Pos2DFinTrac
%     X0: Spot 1 X-position, pixels
%     X1: Spot 2 X-position, pixels
%     Y0: Spot 1 Y-position, pixels
%     Y1: Spot 2 Y-position, pixels
%     Bck: Background level (peak of wide Gaussian 
%     Pk1: Spot 1 Peak, a.u.
%     Pk2: Spot 2 Peak, a.u.
%     spots: spot number present, according to rejection criteria
%-------------------------------------------------------------------------

    r=kymoprops.width;
    LL=kymoprops.duration;
     [~,lp]=size(ReplicationCluster);  
     
     disp('Clusters left = ')
     for i=1:lp;   % %for each cluster
         fprintf(num2str(lp-i,'%03.0f'))
         
         ThisCluster=ReplicationCluster(i);
         startfr=ceil(ThisCluster.PosKyTracCom.frames(1)); %go to start point
         stopfr=ceil(ThisCluster(1,end).PosKyTracCom.frames(end));

         lfr=length(ThisCluster.PosKyTracCom.frames);

         composses=zeros(lfr,1);
         spot1posses=zeros(lfr,1);
         spot2posses=zeros(lfr,1);
         
         spot1amplis=zeros(lfr,1);
         spot1sigmas=zeros(lfr,1);
         
         spot2amplis=zeros(lfr,1);
         spot2sigmas=zeros(lfr,1);
         
         backgrounds=zeros(lfr,1);

         lastpos=ThisCluster.PosKyTracCom.clickpos(1);  %starting position based on manual clicking

         for fr=1:lfr  %for all frames where spot exists:
             
             %We make an estimate based on the overall observed doubling time 
             %measured a priori, for example by checking the kymographs in
             %image J). It is assumed that the channel is filled with
             %growing bacteria up to the last position:
             estpos=lastpos*2^(1/initval.estimateddoublingtime);             
             %estpos=lastpos;
             
             thisfr=ceil(ThisCluster.PosKyTracCom.frames(fr));
             thisline=kymo_FL(thisfr,:);                    %next line (unless last of kymograph)       
             lo=ceil(max(estpos-6, 1)); %select  lo and hi coords
             hi=ceil(min(estpos+6, r)); 
             
             %pre-cook fit
             OriSpotSoi=thisline(lo:hi); 
             MxS=max(OriSpotSoi);
             MnS=min(OriSpotSoi);
             
             SpotSoi=(OriSpotSoi-MnS)/(MxS-MnS); %normalized, background subtracted
           
             %Do Center of mass on normalized line; masked by Gaussian----
             SpotSoi=GaussMask(SpotSoi,0.5);
             [~,comc,~]=Get_1DCOM(SpotSoi); 
             composses(fr)=estpos+comc;   
             lastpos=estpos+comc;       %update last known position
         end
              ReplicationCluster(i).PosKyTracCom.trackpos=composses;
              fprintf('\b\b\b')
     end
     