 function [Ecoli,bacprops]=Division_CheckFluo(Ecoli,framestate,bacprops,bacno,kymo);
     %check if spots are visible in fluorescence channel
    Rightposses=Ecoli(bacno).right;
    lp=length(Rightposses);
    
    
    div=0;
    switch bacprops.sim
        case 1
            spot=0;              %divion condition, simple sim
        case 0
            bacfrno=find(framestate.livingbac==bacno);       %this number gives the index to call this bacterium
            
            %Septum detection %%%%%%%%%%%%%%%%%%%% 
             
            Leftposses=Ecoli(bacno).left;
            lp=length(Leftposses); 
            lastpos_L=Leftposses(lp);
            

            Rightposses=Ecoli(bacno).right;
            lp=length(Rightposses); 
            lastpos_R=Rightposses(lp);
            

            [r,LL]=size(kymo);
            thisfr=min(r,framestate.frame);
            thisline=kymo(thisfr,:);                    %next line (unless last of kymograph)
%             figure;
%             plot(thisline)
%             [~]=ginput(1);
%             close(gcf);
            lo=ceil(max(lastpos_L, 1)); lo=min(LL-1,lo);
            hi=ceil(min(lastpos_R, LL)); hi=max(lo+1, hi);
            
            Fluosoi=thisline(lo:hi);
            Fluosoi=Fluosoi-min(Fluosoi);
            %Fluosoi=GaussMask(Fluosoi,0.5);
            
            tr=8*bacprops.FLtreshold;
            spot=Ifspot(Fluosoi,tr);
            %div=0;
%             nextfr;
%             div;
    end
    FL=Ecoli(bacno).SpotChan1;  
    FL=[FL spot];
    Ecoli(bacno).SpotChan1=FL;  %add to database
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
   
 end