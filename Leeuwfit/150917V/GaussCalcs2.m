%% Load 2D Gaussian fit data and perform calculations.

clear all
close all
clc

% Tus - dif project
expno='001_DnaN_TUS_dif_30122014_difsignal';
%initval=A001_Images_Set_Experiment(expno);      
initval.basepath='/Users/rleeuw/Work/Data/141230_dnaN_dif_tus/dnaN_dif_tus_40msExpTime_5minAcqTimeYFP_30msExpCFP_002_C1/';

% oriZ - Dif project
%initval.basepath='/Users/rleeuw/Work/Data/OriZ-Dif_Results/';
%% Define variables

Ncells=21;

T=cell(Ncells+1,1);
d=cell(Ncells+1,1);
S=cell(Ncells,1);
Sd=cell(Ncells,1);
framesd=zeros(Ncells,1);
framesT=zeros(Ncells,1);

BacLife=zeros(Ncells,1); BacLifed=zeros(Ncells,1);
MainPathTus=strcat(initval.basepath,'Stacks/Tus/DataMULTI/');
MainPathdif=strcat(initval.basepath,'Stacks/dif/DataMULTI/');


%% Load data
tic

for i=1:Ncells;
    T{i}=load(strcat(MainPathTus,num2str(i),'.mat'));
    framesT(i)=length(T{i}.x{1}(:,1));
    d{i}=load(strcat(MainPathdif,num2str(i),'.mat'));
    framesd(i)=length(d{i}.x{1}(:,1));
end

Nspots=T{i}.Nspots;
Npointsd=sum(framesd); 
NpointsT=sum(framesT);
%% Calculations and Filtering

MeanIntT=cell(Nspots,1);MeanIntFCT=cell(Nspots,1);StdIntT=cell(Nspots,1);
StdFCIntT=cell(Nspots,1);MeanIntStrT=cell(Nspots,1);StdIntStrT=cell(Nspots,1);
StdIntd=cell(Nspots,1);StdFCIntd=cell(Nspots,1);MeanIntStrd=cell(Nspots,1);
StdIntStrd=cell(Nspots,1);

n=Ncells+1;


%% Some Calculations


for i=1:Ncells
    BacLife(i)=length(T{i}.x{1}(:,1));
    BacLifed(i)=length(d{i}.x{1}(:,1));
end

Ki=cell(Nspots,1);
KiFC=cell(Nspots,1);
Kx=cell(Nspots,1);
Ky=cell(Nspots,1);

Kdi=cell(Nspots,1);
KdiFC=cell(Nspots,1);
Kdx=cell(Nspots,1);
Kdy=cell(Nspots,1);

MeanBacLifeT=round(mean(BacLife));
MeanBacLifed=round(mean(BacLifed));


for i=1:Ncells
    for j=1:Nspots
    S{i}.x{j}(:,1)=imresize(T{i}.x{j}(:,8),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,7)=imresize(T{i}.x{j}(:,7),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,2)=imresize(T{i}.XNorm{j}(:,2),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,4)=imresize(T{i}.XNorm{j}(:,4),[MeanBacLifeT 1],'bilinear');
    
    Sd{i}.x{j}(:,1)=imresize(d{i}.x{j}(:,8),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,7)=imresize(d{i}.x{j}(:,7),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,2)=imresize(d{i}.XNorm{j}(:,2),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,4)=imresize(d{i}.XNorm{j}(:,4),[MeanBacLifed 1],'bilinear');
   
    end
end

%% Spot tracking algorithms

% 0. Detect particles. Done

% 1. Link particles between consecutive frames. (Using Cost Matrix Linking)

%Cutoffs:

% Should be estimated as 1.05 x maximal cost of all previous links.

bCostd=ones(Nspots,1)*0.3;
dCostd=ones(Nspots,1)*0.3;

bCostT=ones(Nspots,1)*0.3;
dCostT=ones(Nspots,1)*0.3;

bMatrixd=diag(bCostd);
dMatrixd=diag(dCostd);

bMatrixT=diag(bCostT);
dMatrixT=diag(dCostT);

Ctotald=[];
CtotalT=[];

Clinkd=cell(Ncells,MeanBacLifed);
ClinkT=cell(Ncells,MeanBacLifed);

VectorxTd=cell(Ncells,MeanBacLifed);
VectoryTd=cell(Ncells,MeanBacLifed);

VectorxTT=cell(Ncells,MeanBacLifed);
VectoryTT=cell(Ncells,MeanBacLifed);

VectorXTd=cell(Ncells,MeanBacLifed);
VectorYTd=cell(Ncells,MeanBacLifed);

VectorXTT=cell(Ncells,MeanBacLifed);
VectorYTT=cell(Ncells,MeanBacLifed);

VectorITd=cell(Ncells,MeanBacLifed);
VectorITT=cell(Ncells,MeanBacLifed);

VectoriTd=cell(Ncells,MeanBacLifed);
VectoriTT=cell(Ncells,MeanBacLifed);

AuxiliaryMatrixd=cell(Ncells,MeanBacLifed);
AuxiliaryMatrixT=cell(Ncells,MeanBacLifed);

%define initial costs for first frame

% Think about a way to test: 1. constant costs with certain value 2.
% dependency on previous assignment (e.g. max of all previous) 3.
% alternative costs. 
% FIGURE : Displacement of spots vs. number of frames a particle is
% tracked.

for i=1:Ncells
    for t=1:MeanBacLifed-1
        for j=1:Nspots
            for k=1:Nspots
                
                % Cost for linking particle i in fame t to particle j in
                % frame t+1.       
                
                Clinkd{i,t}(j,k)=(sqrt(Sd{i}.x{j}(t,2).^2+Sd{i}.x{j}(t,4).^2)- ...
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
        
        
                AuxiliaryMatrixd{i,t}=0.001*Clinkd{i,t}'; % given the lower cost of the matrix so that it doesn't influence final solution

                Clinkd{i,t}(j+1:j+Nspots,1:k)=bMatrixd; %lower left matrix
                Clinkd{i,t}(1:j,k+1:k+Nspots)=dMatrixd; %upper right matrix
                Clinkd{i,t}(j+1:j+Nspots,k+1:k+Nspots)=AuxiliaryMatrixd{i,t};
                
                Ctotald=[Ctotald Clinkd{i,t}];
                
                bCostd=ones(Nspots,1)*max([bCostd(1) max(Ctotald(:))]); % cost for allowing particles in frame t+1 to get linked by nothing in frame t.
                dCostd=ones(Nspots,1)*max([dCostd(1) max(Ctotald(:))]); % cost for allowing particles in frame t to link to nothing in frame t+1.
                
                bMatrixd=diag(bCostd);
                dMatrixd=diag(dCostd);

                
                Clinkd{i,t}(Clinkd{i,t}==0)=NaN;
                
                [Amind{i,t},Costd{i,t}]=LAP(Clinkd{i,t}); %Linear Assignment Problem (LAP)
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

for i=1:Ncells
    for t=1:MeanBacLifeT-1
        for j=1:Nspots
            for k=1:Nspots
                
                % Cost for linking particle i in fame t to particle j in
                % frame t+1.       

                
                ClinkT{i,t}(j,k)=(sqrt(S{i}.x{j}(t,2).^2+S{i}.x{j}(t,4).^2)- ...
                    sqrt(S{i}.x{k}(t+1,2).^2+S{i}.x{k}(t+1,4).^2)).^2;
               

                %Clinkd{i,t}(j+Nspots,k+Nspots)=Clinkd{i,t}(k,j); %lower right matrix (transpose of upper left)
                
                % t is the frame number
                % i is the cell number
                % j is the spot number
               
            end
            
                    VectorxTT{i,t}=[VectorxTT{i,t} S{i}.x{j}(t+1,2)]; %Vectors used for matrix product with the LAP solution
                    VectoryTT{i,t}=[VectoryTT{i,t} S{i}.x{j}(t+1,4)]; %This will transfer the position data to according spots
                    VectoriTT{i,t}=[VectoriTT{i,t} S{i}.x{j}(t+1,1)];
                    
        end
        
        
                AuxiliaryMatrixT{i,t}=0.001*ClinkT{i,t}'; % given the lower cost of the matrix so that it doesn't influence final solution
                
                ClinkT{i,t}(j+1:j+Nspots,1:k)=bMatrixT; %lower left matrix
                ClinkT{i,t}(1:j,k+1:k+Nspots)=dMatrixT; %upper right matrix
                ClinkT{i,t}(j+1:j+Nspots,k+1:k+Nspots)=AuxiliaryMatrixd{i,t};
                
                CtotalT=[CtotalT ClinkT{i,t}];
                               
                bCostT=ones(Nspots,1)*max([bCostT(1) max(CtotalT(:))]); % cost for allowing particles in frame t+1 to get linked by nothing in frame t.
                dCostT=ones(Nspots,1)*max([dCostT(1) max(CtotalT(:))]); % cost for allowing particles in frame t to link to nothing in frame t+1.
                
                bMatrixT=diag(bCostT);
                dMatrixT=diag(dCostT);
                
                ClinkT{i,t}(ClinkT{i,t}==0)=NaN;
                
                AminT{i,t}=LAP(ClinkT{i,t}); %Linear Assignment Problem (LAP)
        
                %Reassign x-position correspondingly (and corresp
                %intensities)
                %First 'double' the size
                                   
                VectorxTT{i,t}=[VectorxTT{i,t} VectorxTT{i,t}];
                VectoryTT{i,t}=[VectoryTT{i,t} VectoryTT{i,t}];
                VectoriTT{i,t}=[VectoriTT{i,t} VectoriTT{i,t}];
                
                %Matrix Product
               
                VectorXTT{i,t}=AminT{i,t}*VectorxTT{i,t}';
                VectorYTT{i,t}=AminT{i,t}*VectoryTT{i,t}';
                VectorITT{i,t}=AminT{i,t}*VectoriTT{i,t}';
                
                %Correspond results to spots
                
                for j=1:Nspots
                                        
                    S{i}.x{j}(t+1,2)=VectorXTT{i,t}(j);
                    S{i}.x{j}(t+1,4)=VectorYTT{i,t}(j);
                    S{i}.x{j}(t+1,1)=VectorITT{i,t}(j);
                                       
                end
    end
end

% To do 2. close gaps and capture merging and splitting events. (Using Cost Matrix Gap Closing, merging, splitting.)


%% Construct K for taking means of elements at same time point

for i=1:Ncells
    for j=1:Nspots
        
        Ki{j}(:,i)=S{i}.x{j}(:,1);
        KiFC{j}(:,i)=S{i}.x{j}(:,7);

        Kx{j}(:,i)=S{i}.x{j}(:,2);
        Ky{j}(:,i)=S{i}.x{j}(:,4);
        
    end
end



for i=1:Ncells
    for j=1:Nspots
        
        Kdi{j}(:,i)=Sd{i}.x{j}(:,1);
        KdiFC{j}(:,i)=Sd{i}.x{j}(:,7);
        
        Kdx{j}(:,i)=Sd{i}.x{j}(:,2);
        Kdy{j}(:,i)=Sd{i}.x{j}(:,4);
        
    end
end

%% Spot position and intensity filtering. Combination of Spots. --> Within K.
% could be combined within the cost function criteria?

DeltaXcost=0.1^2; % spot cost criterium (10% of image) 
                    % (As a function of corresponding image size later)

DummyVector=ones(Nspots,1)';
VectorX=[];
VectorI=[];
        
for t=1:MeanBacLifed-1
    for i=1:Ncells
        
        for j=1:Nspots; 
            
                    VectorX=[VectorX Sd{i}.x{j}(t,2)];
                    VectorI=[VectorI Sd{i}.x{j}(t,1)]; %These are used to combine later
                    
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
    
    DummyX=VectorX;
    DummyI=VectorI;
    
    if ~isempty(I) 
    
    for l=1:length(I)
        for n=I(l);
            for k=J(l);
            %COM position with intensity weight
            VectorX(n)=(VectorI(n)*DummyX(n)+VectorI(k)*DummyX(k))/(VectorI(n)+VectorI(k));
            VectorX(k)=(VectorI(n)*DummyX(n)+VectorI(k)*DummyX(k))/(VectorI(n)+VectorI(k));
            %Intensities are summed of combined spots
            VectorI(n)=(DummyI(n)+DummyI(k));
            VectorI(k)=(DummyI(n)+DummyI(k));
            end
        end
    end
    end
    
    for j=1:Nspots %This loop can be improved for speed
    Sd{i}.x{j}(t,2)=VectorX(j);
    Sd{i}.x{j}(t,1)=VectorI(j);
    end
    
    VectorX=[];
    VectorI=[];
    clear I
    clear J
    end
end



%% 
M=cell(Nspots,1); Mratio=zeros(Nspots,1); Mratiostd=zeros(Nspots,1);
Md=cell(Nspots,1); Mdratio=zeros(Nspots,1); Mdratiostd=zeros(Nspots,1);

for j=1:Nspots
        
        M{j}(:,1)=mean(Ki{j},2); % mean integrated intensity
        M{j}(:,2)=mean(Kx{j},2); % mean x position
        M{j}(:,3)=std(Kx{j},1,2); % std x position
        M{j}(:,4)=mean(Ky{j},2); % mean y position
        M{j}(:,5)=std(Ky{j},1,2); % std y position
        M{j}(:,6)=std(Ki{j},1,2); % std intensities
        M{j}(:,7)=mean(KiFC{j},2); % mean full cell integrated intensity
        M{j}(:,8)=std(KiFC{j},1,2); % std full cell integrated intensity
        
        Md{j}(:,1)=mean(Kdi{j},2);
        Md{j}(:,2)=mean(Kdx{j},2);
        Md{j}(:,3)=std(Kdx{j},1,2);
        Md{j}(:,4)=mean(Kdy{j},2);
        Md{j}(:,5)=std(Kdy{j},1,2);
        Md{j}(:,6)=std(Kdi{j},1,2);
        Md{j}(:,7)=mean(KdiFC{j},2);
        Md{j}(:,8)=std(KdiFC{j},1,2);

end

%% Filtering

IupboundT=200000;
Iupboundd=150000;

n=Ncells+1;

for i=1:Ncells
    for j=1:T{i}.Nspots
    
% Filter intensities         

    S{i}.x{j}(S{i}.x{j}>=IupboundT)=0;
    Sd{i}.x{j}(Sd{i}.x{j}>=Iupboundd)=0;
    
% Initiate Spot Integrated Intensity Cells

T{n}.IntI{j}=S{1}.x{j}(:,1);
d{n}.IntI{j}=Sd{1}.x{j}(:,1);

% Initiate Full Cell Intensity Cells

T{n}.FCI{j}=S{1}.x{j}(:,7);
d{n}.FCI{j}=Sd{1}.x{j}(:,7);

% Initiate normalised X Position Cells

T{n}.X{j}=S{1}.x{j}(:,2);
d{n}.X{j}=Sd{1}.x{j}(:,2);

% Initiate normalised Y Position Cells

T{n}.Y{j}=S{1}.x{j}(:,4); 
d{n}.Y{j}=Sd{1}.x{j}(:,4);

   end
end

% concatenate matrices to form 

for i=2:Ncells
    for j=1:Nspots

T{n}.FCI{j}=cat(1,T{n}.FCI{j},S{i}.x{j}(:,7));
d{n}.FCI{j}=cat(1,d{n}.FCI{j},Sd{i}.x{j}(:,7));

T{n}.X{j}=cat(1,T{n}.X{j},S{i}.x{j}(:,2));
d{n}.X{j}=cat(1,d{n}.X{j},Sd{i}.x{j}(:,2));

T{n}.Y{j}=cat(1,T{n}.Y{j},S{i}.x{j}(:,4));
d{n}.Y{j}=cat(1,d{n}.Y{j},Sd{i}.x{j}(:,4));

T{n}.IntI{j}=cat(1,T{n}.IntI{j},S{i}.x{j}(:,1));
d{n}.IntI{j}=cat(1,d{n}.IntI{j},Sd{i}.x{j}(:,1));
    end
end

TotCellsStr=sprintf('Ncells = %d',Ncells);

for j=1:Nspots
    
MeanIntT{j}=mean(T{n}.IntI{j});
MeanIntFCT{j}=mean(T{n}.FCI{j});
MeanIntd{j}=mean(d{n}.IntI{j});
MeanIntFCd{j}=mean(d{n}.FCI{j});

StdIntT{j}=std(T{n}.IntI{j});
StdFCIntT{j}=std(T{n}.FCI{j});

MeanIntStrT{j}=sprintf('Mean = %g',MeanIntT{j});
StdIntStrT{j}=sprintf('std = %g',StdIntT{j});

StdIntd{j}=std(d{n}.IntI{j});
StdFCIntd{j}=std(d{n}.FCI{j});

MeanIntStrd{j}=sprintf('Mean = %g',MeanIntd{j});
StdIntStrd{j}=sprintf('std = %g',StdIntd{j});

end

%% Importing the real lifetime data from 'original' Bacpics
MainPathTusOri=strcat(initval.basepath,'Stacks/Tus/DataMULTI/');
MainPathdifOri=strcat(initval.basepath,'Stacks/dif/DataMULTI/');

for i=1:Ncells;
    Tori{i}=load(strcat(MainPathTusOri,num2str(i),'.mat'));
    dori{i}=load(strcat(MainPathdifOri,num2str(i),'.mat'));
    CelllifeT(i)=length(Tori{i}.XNorm(:,1));
    Celllifed(i)=length(dori{i}.XNorm(:,1));
end

MeanCellLifed=mean(Celllifed);
MeanCellLifeT=mean(CelllifeT);

toc
