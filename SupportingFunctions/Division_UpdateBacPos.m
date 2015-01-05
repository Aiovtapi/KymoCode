function Ecoli=Division_UpdateBacPos(Ecoli,framestate,bacprops,bacno,kymo);
%This function updates&elongates the edge position data of a bacterium in a channel
    %Leftposses=Bac{4};
    
    
%     Division = repmat(struct(...
%             'name',0,...            single number indicating family tree position
%             'frames',[1],...         frames            
%             'left',[],...           array containing left edge coordinates from birth until split
%             'right',[],...          array containing right edge coordinates from birth until split
%             'fate','nonexistent',... %'alive' 'exit'  'divided' or 'endofMeas'   
%             'firstfr',[1],...           birth
%             'lastfr',[1],...  %division/end
%             'firstleft',[1],...           birth
%             'firstright',[1],...           birth
%             'lastleft',[1],...  %division/end
%             'lastright',[1]),n_tracks,1);  %division/end

         %Note: births for this set are put at zero, by definition
    bacfrno=find(framestate.livingbac==bacno);       %this number gives the index to call this bacterium
    leftmost=min(framestate.leftposses);   %this is the first position of the first bacterium
    order=framestate.order(bacfrno)-1;      %this number tells how many  bacteria are to the left
    frame=framestate.frame;                                         %(used for predicting the growth)
    
    ThisBac=Ecoli(bacno);
    
    Frames=ThisBac.frames;
    
     
    Leftposses=ThisBac.left; 
    prednextpos_L=Leftposses(end)+bacprops.SimBacGro*order;       %predicted next point
    

    Rightposses=ThisBac.right;
    prednextpos_R=Rightposses(end)+bacprops.SimBacGro*(order+1);     %predicted next point

    nextpos_R=prednextpos_R;
    nextpos_L=prednextpos_L;
    
    if bacprops.sim==0                           %add real data input
        nextline=kymo(framestate.frame,:);
%         figure;
%         plot(nextline)
%         [~]=ginput(1);
%         close(gcf);
        
        LL=length(nextline);
        %get line piece around peak, left position
        lo=ceil(max(prednextpos_L-6, 1)); lo=min(LL-1,lo);
        hi=ceil(min(prednextpos_L+6, LL)); hi=max(lo+1, hi);
        spotsoi_L=nextline(lo:hi);
        spotsoi_L=GaussMask(spotsoi_L,1);
        
        [~,comc_L,~]=Get_1DCOM(spotsoi_L); %Get_1DCom;
        
        %get line piece around peak, right position
        lo=ceil(max(prednextpos_R-6, 1)); lo=min(LL-1,lo);
        hi=ceil(min(prednextpos_R+6, LL)); hi=max(lo+1, hi);
        spotsoi_R=nextline(lo:hi);
        spotsoi_R=GaussMask(spotsoi_R,1);
        [~,comc_R,~]=Get_1DCOM(spotsoi_R); %Get_1DCom;


        nextpos_R=prednextpos_R+comc_R;
        nextpos_L=prednextpos_L+comc_L;
    end
    
    
    Frames=[Frames ; frame];
    Ecoli(bacno).frames=Frames;
    Rightposses=[Rightposses ; nextpos_R];  %update pos
    Ecoli(bacno).right=Rightposses;   
    Leftposses=[Leftposses ; nextpos_L];  %update pos
    Ecoli(bacno).left=Leftposses;
end          


