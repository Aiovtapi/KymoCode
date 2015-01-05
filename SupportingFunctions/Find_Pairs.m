
function pc2=Find_Pairs(pc,bb,initval)
%This function find for each position nearby candidates . 
%JacobKers12
    pc2=0*pc;
    linegrid=10;          %array to be used for interpolation
    pad=3;                  %same, to get 'over the top'
    [lp,~]=size(pc);  c=0;
    for i=1:lp
        x0=pc(i,1); y0=pc(i,2); I=pc(i,3);    %xyI coordinates
        dist=((pc(:,1)-x0).^2+(pc(:,2)-y0).^2).^0.5;
        sel=find((dist<initval.bacteriallength)&dist~=0);      % possible pairs
        if ~isempty(sel);       %close neighbours present
            ls=length(sel)
            for j=1:ls
                x1=pc(sel(j),1); 
                y1=pc(sel(j),2);
                
                
                
                xline=(x1-x0+2*pad).*[1:linegrid]/linegrid+x0-pad;
                yline=(y1-y0+2*pad).*[1:linegrid]/linegrid+y0-pad;
                
                connection=interp2(bb,yline,xline,'linear',0);
                figure;
                plot(connection,'-o');
                [~]=ginput(1)
                close(gcf);
                dum=1;
            end
            %build a connection line here
            %linex=
        else
             
        end
    end
    pc2=pc2(1:c,:);
end


 %profile=interp2(plane,x_intpol,y_intpol,'linear',0); use for line