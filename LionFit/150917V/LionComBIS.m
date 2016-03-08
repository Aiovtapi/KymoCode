function [Sd] = LionComBIS(Sd,DeltaXcost,MeanBacLifed,Ilb)
%LIONCOMBI Summary of this function goes here
%   Detailed explanation goes here

Ncells=size(Sd,2);

VectorXd2=[];
VectorId2=[];
VectorXstdd=[];

for i=1:Ncells
    for t=1:MeanBacLifed;  
        Nspots=size(Sd{i}.x,2);

            for j=1:Nspots; 

                    VectorXstdd=[VectorXstdd Sd{i}.x{j}(t,3)];
                    VectorXd2=[VectorXd2 Sd{i}.x{j}(t,2)];
                    VectorId2=[VectorId2 Sd{i}.x{j}(t,1)];  % This is to test non-COM combination but Brightest Spot's position is taken.     

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
    
    DummyXd=VectorXd2;
    DummyId=VectorId2;
    
    if ~isempty(I) 
    
    for l=1:length(I)
        for n=I(l);
            for k=J(l);
            
            if DummyId(n)>DummyId(k) 
                VectorXd2(n)=DummyXd(n);
                VectorXd2(k)=0;
                VectorId2(n)=(DummyId(n)+DummyId(k));
                VectorId2(k)=0;
            else
                VectorXd2(n)=0;
                VectorXd2(k)=DummyXd(k);
                VectorId2(k)=(DummyId(n)+DummyId(k));
                VectorId2(n)=0;
            end         
            
            end
        end
    end
    end
    
    % Main shift in position and intensities
    
    for j=1:Nspots %This loop can be improved for speed
       
    Sd{i}.x{j}(t,2)=VectorXd2(j); 
    Sd{i}.x{j}(t,1)=VectorId2(j);
                
        if Sd{i}.x{j}(t,1)<Ilb

            Sd{i}.x{j}(t,1:7)=0;
        end

    end
    
    VectorXd2=[];
    VectorId2=[];
    VectorXstdd=[];
    
    clear I
    clear J
    end
end


end

