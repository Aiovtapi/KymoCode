function [DivClicks,Division]=Division_Cluster_Cleaning(DivClicks,Division);
    %Re-check the database inputs
    %label all bad results as fate='badMeas'
    %cluster should live more than 3 frames
    %firstframe should be smaller than last frame
    %at end, remove 'badMeas'
    %at end, keep only fate: 'dissasembled'

% Division = repmat(struct(...
%             'name',0,...            single ['1010111'] style numer indicating family tree position
%             'frames',[1],...         frames            
%             'left',[],...           array containing left edge coordinates from birth until split
%             'right',[],...          array containing right edge coordinates from birth until split
%             'fate','nonexistent'),n_tracks,1);   %'alive' 'exit'  'divided' or 'endofMeas'   
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bax = struct2cell(DivClicks);
fate=squeeze(Bax(2,:,:));   
subset=find(strcmp(fate, 'nonexistent')~=1);  %keep indices of present bacteria
DivClicks=DivClicks(subset);
Division=Division(subset);
    
    
%Keep complete events%%%%%%%%%%%%%%%%5    
    
Buf = struct2cell(DivClicks);
fate=squeeze(Buf(2,:,:));   
subset=find(strcmp(fate, 'divided')==1);  %keep indices of present bacteria
DivClicks=DivClicks(subset);
Division=Division(subset);


  %Identify bad measurements%%%%%%%%%%%%%%%%%%%  
    [~,LD]=size(DivClicks);
    for i=1:LD
        ThisBac=Division(i);
        frs=ThisBac.PosKyTrac.frames;  lfr=length(frs);
        if lfr<4, DivClicks(i).fate='BadMeas'; end
    end
    