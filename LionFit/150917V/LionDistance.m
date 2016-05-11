function [d,dweighted] = LionDistance(S,Sd,LifeT,Lifed)
%LIONDISTANCE Summary of this function goes here
%   Detailed explanation goes here

Ncells=size(Sd,1);

for i=1:Ncells
        LifeT=size(S{i}.x{1},1);
        Lifed=size(Sd{i}.x{1},1);
        NspotsD=size(Sd{i}.x,2);
        NspotsT=size(S{i}.x,2);
        for j=1:NspotsD;        
            %rows: spot position in row number frame.
            Dx{i}(:,j)=Sd{i}.x{j}(:,2);
            DI{i}(:,j)=Sd{i}.x{j}(:,1);  
        end
        
        for j=1:NspotsT
        Tx{i}(:,j)=S{i}.x{j}(:,2);
        TI{i}(:,j)=S{i}.x{j}(:,1);
        end

        %Mean of x of spot in per frame
        TxMeanSpots{i}=nanmean(Tx{i},2);
        DxMeanSpots{i}=nanmean(Dx{i},2);

        TxMeanSpotsW{i}=nansum(TI{i}.*Tx{i},2)./(nansum(TI{i},2));
        DxMeanSpotsW{i}=nansum(DI{i}.*Dx{i},2)./(nansum(DI{i},2));
        
        if Lifed>LifeT
            
            LengthLong=size(DxMeanSpots{i},1);
            LengthShort=size(TxMeanSpots{i},1);
            
            for n=1:LengthShort
                
                Ln(n)=round(n*(LengthLong/LengthShort));
                L=length(DxMeanSpots{i});
                
                if Ln(n)+1<L && Ln(n)-1>0
                DxMeanSpotsShortened{i}(n,1)=mean(DxMeanSpots{i}(Ln(n)-1:1:Ln(n)+1));
                DxMeanSpotsWShortened{i}(n,1)=mean(DxMeanSpotsW{i}(Ln(n)-1:1:Ln(n)+1));
                elseif Ln(n)+1>L
                DxMeanSpotsShortened{i}(n,1)=mean(DxMeanSpots{i}(Ln(n)-2:1:L));
                DxMeanSpotsWShortened{i}(n,1)=mean(DxMeanSpotsW{i}(Ln(n)-2:1:L));
                elseif Ln(n)-1<0
                DxMeanSpotsShortened{i}(n,1)=mean(DxMeanSpots{i}(1:1:Ln(n)+2));
                DxMeanSpotsWShortened{i}(n,1)=mean(DxMeanSpotsW{i}(1:1:Ln(n)+2));
                end
            end 
            d{i}=abs(DxMeanSpotsShortened{i}-TxMeanSpots{i});
            dweighted{i}=abs(DxMeanSpotsWShortened{i}-TxMeanSpotsW{i});
        else
            LengthLong=size(TxMeanSpots{i},1);
            LengthShort=size(DxMeanSpots{i},1);
            
            for n=1:LengthShort
                
                Ln(n)=round(n*(LengthLong/LengthShort));
                L=length(TxMeanSpots{i});
                
                if Ln(n)+1<L && Ln(n)-1>0
                TxMeanSpotsShortened{i}(n,1)=mean(TxMeanSpots{i}(Ln(n)-1:1:Ln(n)+1));
                TxMeanSpotsWShortened{i}(n,1)=mean(TxMeanSpotsW{i}(Ln(n)-1:1:Ln(n)+1));
                elseif Ln(n)+1>L
                TxMeanSpotsShortened{i}(n,1)=mean(TxMeanSpots{i}(Ln(n)-2:1:L));
                TxMeanSpotsWShortened{i}(n,1)=mean(TxMeanSpotsW{i}(Ln(n)-2:1:L));
                elseif Ln(n)-1<0
                TxMeanSpotsShortened{i}(n,1)=mean(TxMeanSpots{i}(1:1:Ln(n)+2));
                TxMeanSpotsWShortened{i}(n,1)=mean(TxMeanSpotsW{i}(1:1:Ln(n)+2));
                end
                
            end
            
            d{i}=abs(TxMeanSpotsShortened{i}-DxMeanSpots{i});
            dweighted{i}=abs(TxMeanSpotsWShortened{i}-DxMeanSpotsW{i});
        end
        
end

end

