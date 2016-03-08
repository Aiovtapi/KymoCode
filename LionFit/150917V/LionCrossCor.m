function [r] = LionCrossCor(S,Sd,MeanBacLifed,MeanBacLifeT)
%LIONCROSSCOR Summary of this function goes here
%   Detailed explanation goes here

Ncells=size(Sd,1);
Nspots=size(Sd{1}.x,2);

    
%cross-correlation is function of time, x and y.
%(but we don't have equal number of points)
%This is a big problem.

%I will go through the channel with least frames.
%Then find spots in nearest frame in the other channel.


for i=1:Ncells;
        for j=1:Nspots;
            
            %rows: spot position in row number frame.
            Tx{i}(:,j)=S{i}.x{j}(:,2);
            Dx{i}(:,j)=Sd{i}.x{j}(:,2);
            
        end
    
        %Mean of x of spot in per frame
        TxMeanSpots{i}=nanmean(Tx{i},2);
        DxMeanSpots{i}=nanmean(Dx{i},2);
    
        TxMeanOverall{i}=nanmean(TxMeanSpots{i},1);
        DxMeanOverall{i}=nanmean(DxMeanSpots{i},1);
        
        if MeanBacLifed>MeanBacLifeT
            
            LengthLong=size(DxMeanSpots{i},1);
            LengthShort=size(TxMeanSpots{i},1);
            

                    for n=1:LengthShort
                        Ln(n)=round(n*(LengthLong/LengthShort));
                        DxMeanSpotsShortened{i}(n)=DxMeanSpots{i}(Ln);
                    end
                    
                    r{i}=sum(TxMeanSpots{i}*DxMeanSpots{i});
        else
            LengthLong=size(DxMeanSpots{i},1);
            LengthShort=size(TxMeanSpots{i},1);            
        end
        

end








end

