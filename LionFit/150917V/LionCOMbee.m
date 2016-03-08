function [Sd,Intensities] = LionCOMbee(Sd,DeltaXcost,MeanBacLifed,Ilb)
%LIONCOMBI Summary of this function goes here
%   Detailed explanation goes here

Ncells=size(Sd,1);
Nspots=size(Sd{1}.x,2);

VectorXd=[];
VectorId=[];
VectorXstdd=[];


        
for t=1:MeanBacLifed-1
    for i=1:Ncells
        
        for j=1:Nspots; 
            
                    VectorXd=[VectorXd Sd{i}.x{j}(t,2)];
                    VectorId=[VectorId Sd{i}.x{j}(t,1)]; %These are used to combine later
                    VectorXstdd=[VectorXstdd Sd{i}.x{j}(t,3)];
                

            for k=1:Nspots;  
                
                %New combination matrix
                Ccombd{i,t}(j,k)=(sqrt(Sd{i}.x{j}(t,2).^2+Sd{i}.x{j}(t,4).^2)- ...
                    sqrt(Sd{i}.x{k}(t,2).^2+Sd{i}.x{k}(t,4).^2)).^2;

            end
        end
                        
    Closed{i,t}=Ccombd{i,t}<DeltaXcost;
        
    Closed{i,t}=Closed{i,t}-diag(ones(Nspots,1)); %remove diagonals, they are always one.
        
    [I,J]=find(tril(Closed{i,t})==1); 
    %I and J are indices of the nonzero non-diagonal elements
    %which indicate a combination.
    
    DummyXd=VectorXd;
    DummyId=VectorId;
    
    if ~isempty(I) 
    
    for l=1:length(I)
        for n=I(l);
            for k=J(l);

            %COM position with intensity weight
            VectorXd(n)=(DummyId(n)*DummyXd(n)+DummyId(k)*DummyXd(k))/(DummyId(n)+DummyId(k));
            VectorXd(k)=(DummyId(n)*DummyXd(n)+DummyId(k)*DummyXd(k))/(DummyId(n)+DummyId(k));        

            %Intensities are summed of combined spots
            VectorId(n)=(DummyId(n)+DummyId(k));
            VectorId(k)=(DummyId(n)+DummyId(k));

            if DummyId(n)>DummyId(k) % This is for making accurate intensity calculations (avoid double counting)
            VectorId2(n)=(DummyId(n)+DummyId(k));
            VectorId2(k)=NaN;  
            else
            VectorId2(k)=(DummyId(n)+DummyId(k));
            VectorId2(n)=NaN;   
            end
            
            end
        end
    end
    
    for j=1:Nspots %This loop can be improved for speed
    Sd{i}.x{j}(t,2)=VectorXd(j);
    Sd{i}.x{j}(t,1)=VectorId(j);
    
    Intensities{i}.x{j}(t,1)=VectorId2(j);
    
        if Sd{i}.x{j}(t,1)<Ilb
            Sd{i}.x{j}(t,1)=NaN;
            Sd{i}.x{j}(t,2)=NaN;
            
            Intensities{i}.x{j}(t,1)=NaN;
        end
        
    end
    
    VectorXd=[];
    VectorId=[];
    VectorXstdd=[];

    clear I
    clear J
    end
end


end

