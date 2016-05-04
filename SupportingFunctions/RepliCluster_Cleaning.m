function [ReplicationCluster,RepClicks]=RepliCluster_Cleaning(ReplicationCluster,RepClicks)
    %Re-check the database inputs
    %label all bad results as fate='badMeas'
    %cluster should live more than 3 frames
    %firstframe should be smaller than last frame
    %at end, remove 'badMeas'
    %at end, keep only fate: 'dissasembled'
    
    
  %Identify bad measurements%%%%%%%%%%%%%%%%%%%  
    [~,LR]=size(RepClicks);
    for i=1:LR
        ThisRep=RepClicks(i);
        fr1=ThisRep.PosClick.firstframe; fr2=ThisRep.PosClick.lastframe;
        %firstframe should be smaller than last frame
        if fr2-fr1<4, RepClicks(i).fate='BadMeas'; end
        %cluster should live more than 3 frames
        frames=ReplicationCluster(i).PosKyTracCom.frames; lfr=length(frames);
        if lfr<4, RepClicks(i).fate='BadMeas';     end
    end
 
%Keep complete events%%%%%%%%%%%%%%%%5 
if 1
Buf = struct2cell(RepClicks);
fate=squeeze(Buf(2,:,:));   
subset=find(strcmp(fate, 'BadMeas')~=1);  %keep indices of present bacteria
ReplicationCluster=ReplicationCluster(subset);
RepClicks=RepClicks(subset);
end
    