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

Ncells=3;

T=cell(Ncells+1,1);
d=cell(Ncells+1,1);
S=cell(Ncells,1);
Sd=cell(Ncells,1);
framesd=zeros(Ncells,1);
framesT=zeros(Ncells,1);

BacLife=zeros(Ncells,1); BacLifed=zeros(Ncells,1);
MainPathTus=strcat(initval.basepath,'StacksLong/Tus/DataMULTI/');
MainPathdif=strcat(initval.basepath,'StacksLong/dif/DataMULTI/');


%% Load data
tic

for i=1:Ncells;
    T{i}=load(strcat(MainPathTus,num2str(i),'.mat'));
    framesT(i)=length(T{i}.x{1}(:,1));
    d{i}=load(strcat(MainPathdif,num2str(i),'.mat'));
    framesd(i)=length(d{i}.x{1}(:,1));
end

Nspots=T{i}.Nspots;

%% Calculations and Filtering

MeanIntT=cell(Nspots,1);MeanIntFCT=cell(Nspots,1);StdIntT=cell(Nspots,1);
StdFCIntT=cell(Nspots,1);MeanIntStrT=cell(Nspots,1);StdIntStrT=cell(Nspots,1);
StdIntd=cell(Nspots,1);StdFCIntd=cell(Nspots,1);MeanIntStrd=cell(Nspots,1);
StdIntStrd=cell(Nspots,1);

IupboundT=200000;
Iupboundd=150000;

n=Ncells+1;

for i=1:Ncells
    for j=1:T{i}.Nspots
    
% Filter intensities         

    T{i}.x{j}(T{i}.x{j}>=IupboundT)=0;
    d{i}.x{j}(d{i}.x{j}>=Iupboundd)=0;
    
% Initiate Spot Intensity Cells

T{n}.I{j}=T{1}.x{j}(:,1);
d{n}.I{j}=d{1}.x{j}(:,1);

% Initiate Cell Intensity Cells

T{n}.FCI{j}=T{1}.x{j}(:,7);
d{n}.FCI{j}=d{1}.x{j}(:,7);

% Initiate normalised X Position Cells

T{n}.X{j}=T{1}.XNorm{j}(:,2);
d{n}.X{j}=d{1}.XNorm{j}(:,2);

% Initiate normalised Y Position Cells

T{n}.Y{j}=T{1}.XNorm{j}(:,4); 
d{n}.Y{j}=d{1}.XNorm{j}(:,4);

% Initiate Integrated Intensity Cells

T{n}.IntI{j}=T{1}.x{j}(:,8);
d{n}.IntI{j}=d{1}.x{j}(:,8);
   end
end

% concatenate matrices to form 
for i=2:Ncells
    for j=1:Nspots
T{n}.I{j}=cat(1,T{n}.I{j},T{i}.x{j}(:,1));
d{n}.I{j}=cat(1,d{n}.I{j},d{i}.x{j}(:,1));

T{n}.FCI{j}=cat(1,T{n}.FCI{j},T{i}.x{j}(:,7));
d{n}.FCI{j}=cat(1,d{n}.FCI{j},d{i}.x{j}(:,7));

T{n}.X{j}=cat(1,T{n}.X{j},T{i}.XNorm{j}(:,2));
d{n}.X{j}=cat(1,d{n}.X{j},d{i}.XNorm{j}(:,2));

T{n}.Y{j}=cat(1,T{n}.Y{j},T{i}.XNorm{j}(:,4));
d{n}.Y{j}=cat(1,d{n}.Y{j},d{i}.XNorm{j}(:,4));

T{n}.IntI{j}=cat(1,T{n}.IntI{j},T{i}.x{j}(:,8));
d{n}.IntI{j}=cat(1,d{n}.IntI{j},d{i}.x{j}(:,8));
    end
end

%% Some Calculations
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

MeanBacLifeT=mean(BacLife);
MeanBacLifed=mean(BacLifed);


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
    
    radiusd{i}.x{j}=sqrt(Sd{i}.x{j}(:,2).^2+Sd{i}.x{j}(:,4).^2);
    end
end

%% Spot tracking algorithms

% 0. Detect particles. Done

% 1. Link particles between consecutive frames. (Using Cost Matrix Linking)

%Cutoffs:

% Should be estimated as 1.05 x maximal cost of all previous links.

bCost=ones(Nspots,1)*0.1;
dCost=ones(Nspots,1)*0.1;


bMatrix=diag(bCost);
dMatrix=diag(dCost);

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
                
                Clinkd{i,t}(j,k)=(radiusd{i}.x{j}(t)-radiusd{i}.x{k}(t+1)).^2;

                %Clinkd{i,t}(j+Nspots,k+Nspots)=Clinkd{i,t}(k,j); %lower right matrix (transpose of upper left)
                
                % t is the frame number
                % i is the cell number
                % j is the spot number
                
            end
        end
                AuxiliaryMatrix{i,t}=Clinkd{i,t}'; 
                
                Clinkd{i,t}(j+1:j+Nspots,1:k)=bMatrix; %lower left matrix
                Clinkd{i,t}(1:j,k+1:k+Nspots)=dMatrix; %upper right matrix
                Clinkd{i,t}(j+1:j+Nspots,k+1:k+Nspots)=AuxiliaryMatrix{i,t};
                

                Ctotal=Clinkd{i,t};
                
                bCost=ones(Nspots,1)*max([bCost(1) max(Ctotal(:))]); % cost for allowing particles in frame t+1 to get linked by nothing in frame t.
                dCost=ones(Nspots,1)*max([dCost(1) max(Ctotal(:))]); % cost for allowing particles in frame t to link to nothing in frame t+1.
                
                bMatrix=diag(bCost);
                dMatrix=diag(dCost);
                

                
    end
end


% 2. close gaps and capture merging and splitting events. (Using Cost Matrix Gap Closing, merging, splitting.)



%% Construct K for taking means of elements at same time point

for i=1:Ncells
    for j=1:Nspots
        
        Ki{j}(:,i)=S{i}.x{j}(:,1);
        KiFC{j}(:,i)=S{i}.x{j}(:,7);

        Kx{j}(:,i)=S{i}.x{j}(:,2);
        Ky{j}(:,i)=S{i}.x{j}(:,4);
        
    end
end

for n=1:MeanBacLifed
    for i=1:Ncells
        
        Kdi{j}(:,i)=Sd{i}.x{j}(:,1);
        KdiFC{j}(:,i)=Sd{i}.x{j}(:,7);
        Kdx{j}(:,i)=Sd{i}.x{j}(:,2);
        Kdy{j}(:,i)=Sd{i}.x{j}(:,4);
        
    end
end

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
