function [Sd] = LionLink(Sd,MeanBacLifed)
%% Spot tracking algorithms

Ncells=size(Sd,1);
Nspots=size(Sd{1}.x,2);

% 0. Detect particles. Done

% 1. Link particles between consecutive frames. (Using Cost Matrix Linking)

if nargin<2
    %Cutoffs:
Timelength=[];

for i=1:Ncells
    Timelength=[Timelength size(Sd{i}.x{1},1)];
end

MeanBacLifed=max(Timelength);
end

Ctotald=[];

Clinkd=cell(Ncells,MeanBacLifed);

VectorxTd=cell(Ncells,MeanBacLifed);
VectoryTd=cell(Ncells,MeanBacLifed);

VectorXTd=cell(Ncells,MeanBacLifed);
VectorYTd=cell(Ncells,MeanBacLifed);

VectorITd=cell(Ncells,MeanBacLifed);

VectoriTd=cell(Ncells,MeanBacLifed);

AuxiliaryMatrixd=cell(Ncells,MeanBacLifed);

%define initial costs for first frame

% Think about a way to test: 1. constant costs with certain value 2.
% dependency on previous assignment (e.g. max of all previous) 3.
% alternative costs. 
% FIGURE : Displacement of spots vs. number of frames a particle is
% tracked.


for i=1:Ncells
    
    if nargin<2
    TT=size(Sd{i}.x{1},1);
    else
    TT=MeanBacLifed;
    end
    
    for t=1:TT-1
        for j=1:size(Sd{i}.x,2)
            for k=1:size(Sd{i}.x,2)
                
                Nspots=size(Sd{i}.x,2);
                % Cost for linking particle i in frame t to particle j in
                % frame t+1.
                
                Clinkd{i,t}(j,k)=(sqrt(Sd{i}.x{j}(t,2).^2+Sd{i}.x{j}(t,4).^2) - ...
                    sqrt(Sd{i}.x{k}(t+1,2).^2+Sd{i}.x{k}(t+1,4).^2)).^2;

                %Clinkd{i,t}(j+Nspots,k+Nspots)=Clinkd{i,t}(k,j); %lower right matrix (transpose of upper left)
                
                % t is the frame number
                % i is the cell number
                % j is the spot number
            end
                    VectorxTd{i,t}=[VectorxTd{i,t} Sd{i}.x{j}(t+1,2)]; %Vectors used for matrix product with the LAP solution
                    VectoryTd{i,t}=[VectoryTd{i,t} Sd{i}.x{j}(t+1,4)]; %This will transfer the position data to according spots
                    VectoriTd{i,t}=[VectoriTd{i,t} Sd{i}.x{j}(t+1,1)];
                  
        end
        
        % Should be estimated as 1.05 x maximal cost of all previous links.

                bCostd=ones(Nspots,1)*2E5; % This is still to be implemented (hence 20).
                dCostd=ones(Nspots,1)*2E5;

                bMatrixd=diag(bCostd);
                dMatrixd=diag(dCostd);
                
                AuxiliaryMatrixd{i,t}=2E5*Clinkd{i,t}'; % given the lower cost of the matrix so that it doesn't influence final solution

                Clinkd{i,t}(j+1:j+Nspots,1:k)=bMatrixd; %lower left matrix
                Clinkd{i,t}(1:j,k+1:k+Nspots)=dMatrixd; %upper right matrix
                Clinkd{i,t}(j+1:j+Nspots,k+1:k+Nspots)=AuxiliaryMatrixd{i,t};
                
                Ctotald=[Ctotald;Clinkd{i,t}(:)];
                
                bCostd=ones(Nspots,1)*max([bCostd(1) max(Ctotald(:))]); % cost for allowing particles in frame t+1 to get linked by nothing in frame t.
                dCostd=ones(Nspots,1)*max([dCostd(1) max(Ctotald(:))]); % cost for allowing particles in frame t to link to nothing in frame t+1.
                
                bMatrixd=diag(bCostd);
                dMatrixd=diag(dCostd);

                
                Clinkd{i,t}(Clinkd{i,t}==0)=NaN;
                
                %Linear Assignment Problem (LAP)
                
%                 [Amind{i,t},Costd{i,t}]=munkres(Clinkd{i,t}); 
                [Amind{i,t},Costd{i,t}]=LionLAP(Clinkd{i,t});
                
                %Costd is the cost corresponding with the switches. Numbers
                %correspond with the columns.
                
                %Reassign x-position correspondingly (and corresp
                %intensities)
                %First 'double' the size
                    
                VectorxTd{i,t}=[VectorxTd{i,t} VectorxTd{i,t}];
                VectoryTd{i,t}=[VectoryTd{i,t} VectoryTd{i,t}];
                VectoriTd{i,t}=[VectoriTd{i,t} VectoriTd{i,t}];
                
                %Matrix Product
                VectorXTd{i,t}=Amind{i,t}*VectorxTd{i,t}';
                VectorYTd{i,t}=Amind{i,t}*VectoryTd{i,t}';
                VectorITd{i,t}=Amind{i,t}*VectoriTd{i,t}';
                             
                %Correspond results to spots
                
                for j=1:Nspots
                                        
                    Sd{i}.x{j}(t+1,2)=VectorXTd{i,t}(j);
                    Sd{i}.x{j}(t+1,4)=VectorYTd{i,t}(j);
                    Sd{i}.x{j}(t+1,1)=VectorITd{i,t}(j);
                    
                end
                
    end
end

end

