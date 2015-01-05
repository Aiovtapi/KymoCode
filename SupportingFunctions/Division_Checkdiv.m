 function [Ecoli,bacprops]=Division_Checkdiv(Ecoli,framestate,bacprops,bacno,kymo);
     %check if division took place on a specific condition
    
     
     ThisBac=Ecoli(bacno);
    
    Leftposses=ThisBac.left;
    lastpos_L=Leftposses(end);
            
    Rightposses=ThisBac.right;
    lastpos_R=Rightposses(end);
    lp=length(Rightposses);
    div=0;
    switch bacprops.sim
        case 1
            if lastpos_R-lastpos_L>=bacprops.SimMaxL, div=1;, end;           
        case 0
            bacfrno=find(framestate.livingbac==bacno);       %this number gives the index to call this bacterium
            order=framestate.order(bacfrno)-1;      %this number tells how many  bacteria are to the left
             %(used for predicting the growth)
            %Septum detection %%%%%%%%%%%%%%%%%%%% 
             bacgro=0.2;
            
            

            [r,LL]=size(kymo);
            thisfr=min(r,framestate.frame);
            thisline=kymo(thisfr,:);                    %next line (unless last of kymograph)
            lo=ceil(max(lastpos_L, 1)); lo=min(LL-1,lo);
            hi=ceil(min(lastpos_R, LL)); hi=max(lo+1, hi);
            
            Divsoi=thisline(lo:hi);
            Divsoi=Divsoi-min(Divsoi);
            Divsoi=GaussMask(Divsoi,0.5);
            tr=bacprops.treshold;
            div=Ifspot(Divsoi,tr);
            %div=0;
%             nextfr;
%             div;
    end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    if div==1     
        if bacprops.sim==0, 
            [~,comc,~]=Get_1DCOM(Divsoi); 
        else comc=0; %Get_1DCom;
        end
        
        Ecoli(bacno).fate='divided';
        Ecoli(bacno).end=framestate.frame;

        %make a 'left' offspring
        Ecoli(bacprops.count+1).name=double(2*Ecoli(bacno).name);
        Ecoli(bacprops.count+1).frames=framestate.frame+1;
        Ecoli(bacprops.count+1).left=Ecoli(bacno).left(lp);
        Ecoli(bacprops.count+1).right=(Ecoli(bacno).left(lp)+Ecoli(bacno).right(lp))/2+comc;
        Ecoli(bacprops.count+1).fate='alive';

        %make a 'right' offspring
        Ecoli(bacprops.count+2).name=double(2*Ecoli(bacno).name+1);
        Ecoli(bacprops.count+2).frames=framestate.frame+1;
        Ecoli(bacprops.count+2).left=(Ecoli(bacno).left(lp)+Ecoli(bacno).right(lp))/2+comc;
        Ecoli(bacprops.count+2).right=Ecoli(bacno).right(lp);
        Ecoli(bacprops.count+2).fate='alive';

        bacprops.count=bacprops.count+2;

    end
 end