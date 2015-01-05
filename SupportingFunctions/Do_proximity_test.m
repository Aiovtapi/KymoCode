
function pc2=Do_proximity_test(pc,initval)
%This function calculates for each position if a too nearby maximum is
%larger or not. If so, it is discarded. JacobKers11
    pc2=0*pc;
    [lp,~]=size(pc);  c=0;
    for i=1:lp
        x0=pc(i,1); y0=pc(i,2); I=pc(i,3);    %xyI coordinates
        dist=((pc(:,1)-x0).^2+(pc(:,2)-y0).^2).^0.5;
        sel=find((dist<initval.minproximity)&dist~=0);
        if ~isempty(sel);       %close neighbours
            Iprox=pc(sel,3); 
            if I>2*max(Iprox)     %is it -clearly- the local winner?
                c=c+1;
                pc2(c,:)=pc(i,:);
            end
        else
          c=c+1;%isolated maximum
          pc2(c,:)=pc(i,:);    
        end
    end
    pc2=pc2(1:c,:);
end
