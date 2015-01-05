function spots=TimeLine_FindRepClusEvents(ReplicationCluster, ThisBac);
%This function finds which events (in this case, replication clusters) have
%timelines that fall within the range set by time end edge psositions of
%the bacterium

leX_Div=ThisBac.left;
loX_Div=min(leX_Div);


reX_Div=ThisBac.right;
hiX_Div=max(reX_Div);


frs_Div=ThisBac.frames;
loT_Div=min(frs_Div); 
hiT_Div=max(frs_Div);

lc=length(ReplicationCluster);
spots=[];
for i=1:lc
    if abs(log2(ThisBac.name)-log2(ThisBac.name))<2  %only look to direct ancestors, siblings  and successors
    frs_RC=ReplicationCluster(i).frames;
    pos_RC=ReplicationCluster(i).trackpos;
    loT_RC=min(frs_RC); 
    hiT_RC=max(frs_RC);
    loX_RC=min(pos_RC);
    hiX_RC=max(pos_RC);
    
    
    cond1=(hiT_RC>loT_Div);  %spot should end after birth
    cond2=(loT_RC<hiT_Div);  %spot should begin before division
    cond3=(loX_Div<hiX_RC);  %spot should somewhere overlap with Bac
    cond4=(hiX_Div>loX_RC);  %same
    cond=cond1*cond2*cond3*cond4;
    
    inside=0;%final check:the spot should at least somewhere lie in between
    if cond, 
        lfrc=length(frs_RC);       
        for j=1:lfrc
            fr=frs_RC(j);  %current time spot
            x=pos_RC(j);   %current pos spot
            sel=find(ceil(frs_Div)==ceil(fr));  
            left=leX_Div(sel);  %find current edges
            right=reX_Div(sel); %find current edges
            if (x>left&x<right), 
                inside=1;
            end;  %the spot should at least somewhere lie in between
        end
        if inside, spots=[spots i];  end
    end  %these spot traces are overlapping with bact.
  end
end

