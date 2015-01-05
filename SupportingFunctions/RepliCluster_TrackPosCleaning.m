function ReplicationCluster=RepliCluster_TrackPosCleaning(ReplicationCluster);
%Function cleans up tracking results

%ReplicationCluster
% struct array with fields:
%     name
%     frames
%     clickpos
%     trackpos
%     fate
%     spot1pos
%     spot2pos
%     firstpos
%     lastpos
%     firstframe
%     lastframe
%     spot1amp
%     spot1sig
%     spot2amp
%     spot2sig
    
  %Clean up track pos----------------------------------  
    LR=length(ReplicationCluster);
    for j=1:LR
        ThisRep=ReplicationCluster(j);
        ThisRepbuf=ThisRep;
        tr=ThisRep.PosKyTracCom.trackpos; le=length(tr);  
        for i=1:le           
            if ThisRep.PosKyTracGauss.spot2pos(i)<ThisRep.PosKyTracGauss.spot1pos(i), 
                ThisRep.PosKyTracGauss.spot1pos(i)=ThisRepbuf.PosKyTracGauss.spot2pos(i);
                ThisRep.PosKyTracGauss.spot2pos(i)=ThisRepbuf.PosKyTracGauss.spot1pos(i);
                ThisRep.PosKyTracGauss.spot1amp(i)=ThisRepbuf.PosKyTracGauss.spot2amp(i);
                ThisRep.PosKyTracGauss.spot1sig(i)=ThisRepbuf.PosKyTracGauss.spot2sig(i);
                ThisRep.PosKyTracGauss.spot2amp(i)=ThisRepbuf.PosKyTracGauss.spot1amp(i);
                ThisRep.PosKyTracGauss.spot2sig(i)=ThisRepbuf.PosKyTracGauss.spot1sig(i);

            end
        end
        ReplicationCluster(j)=ThisRep;
        %plot(tim, csp2-csp1, '-ro'); hold off;
        %[~]=ginput(1);
    end
 

    