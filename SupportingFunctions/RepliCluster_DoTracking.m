    function ReplicationCluster=RepliCluster_DoTracking(ReplicationCluster,kymoprops,kymo_FL);
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
     for i=1:lp   % %for each cluster
         lp-i
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
             thisfr=ceil(ThisCluster.PosKyTracCom.frames(fr));
             thisline=kymo_FL(thisfr,:);                    %next line (unless last of kymograph)       
             lo=ceil(max(lastpos-6, 1)); %select  lo and hi coords
             hi=ceil(min(lastpos+6, r)); 
             
             %pre-cook fit
             OriSpotSoi=thisline(lo:hi); 
             MxS=max(OriSpotSoi);
             MnS=min(OriSpotSoi);
             
             SpotSoi=(OriSpotSoi-MnS)/(MxS-MnS); %normalized, background subtracted

             %2x 1D GaussFit%%%%%%%%%%%%%%%%%%%%%5
             modus='fixpsf'; %'fixpsf';
             est.psf=1.2;
             c=length(SpotSoi);
             est.N=0.5;  
             est.x0=lo+c/2 ;%-c/4;
             est.x1=lo+c/2; %+c/4;

            %2x 1D GaussFit on %%%%%%%%%%%%%%%%%%%%%5
             plotit=0;
             x=(lo:1:hi);
             switch modus
                 case 'fixpsf'
                    paramsF=MLE_Two1D_Gaussian_RealDataFixPSF(x,SpotSoi,est);
                    %recalculate parameters to: 
                    s1=paramsF(1);  
                    N1=(MxS-MnS)*paramsF(4)*(est.psf*(2*pi)^0.5);  %integrated intensity 
                    sig1=est.psf;
                    s2=paramsF(2);  
                    N2=(MxS-MnS)*paramsF(5)*(est.psf*(2*pi)^0.5); 
                    sig2=est.psf;
                    background=MnS;
                 case 'fitpsf'                
                    paramsF=MLE_Two1D_Gaussian_RealData(x,SpotSoi,est,plotit);  %[ux1 ux2 sigmaValueX1 sigmaValueX2 b N1 N2];
                    s1=paramsF(1);  
                    N1=(MxS-MnS)*paramsF(6)*(est.psf*(2*pi)^0.5); 
                    sig1=paramsF(3);
                    s2=paramsF(2);  
                    N2=(MxS-MnS)*paramsF(7)*(est.psf*(2*pi)^0.5); 
                    sig2=paramsF(4);
                    background=MnS;
             end
            
             spot1posses(fr)=s1;
             if ~isempty(N1), spot1amplis(fr)=N1;, else spot1amplis(fr)=0;, end
             spot1sigmas(fr)=sig1;
             spot2posses(fr)=s2;
             if ~isempty(N2), spot2amplis(fr)=N2;, else spot2amplis(fr)=0;, end
             
             spot2sigmas(fr)=sig2;
             backgrounds(fr)=background;
             
%              sum1=sum(OriSpotSoi) %sum of all pixel values
%              sum2=N1+N2+length(OriSpotSoi)*background
%              dum=1;
             %------------Clean Gauss result from 'unreasonable' values
             check=1;
             if check
                 %Spots in fit range?
                 if s1>lo&s1<hi, 
                     spot1posses(fr)=s1;
                     spot1amplis(fr)=N1;
                     spot1sigmas(fr)=sig1;
                     
                 else
                     spot1posses(fr)=NaN;
                     spot1amplis(fr)=NaN;
                     spot1sigmas(fr)=NaN;
                 end      
                 %%%%%%%%%%%%%%%%%%%
                  if (s2>lo&s2<hi), 
                     spot2posses(fr)=s2;
                     spot2amplis(fr)=N2;
                     spot2sigmas(fr)=sig2;
                  else
                     spot2posses(fr)=NaN;
                     spot2amplis(fr)=NaN;
                     spot2sigmas(fr)=NaN;        
                  end
                 %Spots strong enough?
                 if N1/N2<0.10; 
                     spot1posses(fr)=NaN;
                     spot1amplis(fr)=NaN;
                     spot1sigmas(fr)=NaN;
                 end
                 if N2/N1<0.10; 
                     spot2posses(fr)=NaN;
                     spot2amplis(fr)=NaN;
                     spot2sigmas(fr)=NaN;   
                 end
             else
                 spot1posses(fr)=s1;
                 spot1amplis(fr)=N1;
                 spot1sigmas(fr)=sig1;
                 spot2posses(fr)=s2;
                 spot2amplis(fr)=N2;
                 spot2sigmas(fr)=sig2;
             end
             %--------------------------------------------------
             
             %Do Center of mass on normalized line; masked by Gaussian----
             SpotSoi=GaussMask(SpotSoi,0.5);
            [~,comc,~]=Get_1DCOM(SpotSoi); 
            composses(fr)=lastpos+comc;   
            lastpos=lastpos+comc;       %update last known position
         end
         
         dum=1;

              ReplicationCluster(i).PosKyTracCom.trackpos=composses;
              ReplicationCluster(i).PosKyTracGauss.spot1pos=spot1posses;
              ReplicationCluster(i).PosKyTracGauss.spot1amp=spot1amplis;
              ReplicationCluster(i).PosKyTracGauss.spot1sig=spot1sigmas;
              ReplicationCluster(i).PosKyTracGauss.spot2pos=spot2posses;
              ReplicationCluster(i).PosKyTracGauss.spot2amp=spot2amplis;
              ReplicationCluster(i).PosKyTracGauss.spot2sig=spot2sigmas; 
              ReplicationCluster(i).PosKyTracGauss.background=background;
         dum=1;
     end
     