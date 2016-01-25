%% Load 2D Gaussian fit data and perform calculations.

clear all
close all
clc

% Tus - dif project
 expno='001_DnaN_TUS_dif_30122014_difsignal';
 initval=A001_Images_Set_Experiment(expno); %define your paths and files
 
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

Nspots=T{i}.NSpots;
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
    for j=1:size(d{i}.x,2) 
    Sd{i}.x{j}(:,1)=imresize(d{i}.x{j}(:,8),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,7)=imresize(d{i}.x{j}(:,7),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,3)=imresize(d{i}.x{j}(:,3),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,2)=imresize(d{i}.XNorm{j}(:,2),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,4)=imresize(d{i}.XNorm{j}(:,4),[MeanBacLifed 1],'bilinear');
    end
    
    for j=1:size(T{i}.x,2)
    S{i}.x{j}(:,1)=imresize(T{i}.x{j}(:,8),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,7)=imresize(T{i}.x{j}(:,7),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,3)=imresize(T{i}.x{j}(:,3),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,2)=imresize(T{i}.XNorm{j}(:,2),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,4)=imresize(T{i}.XNorm{j}(:,4),[MeanBacLifeT 1],'bilinear');
    end
end

%% Spot tracking + linking algorithms

%1. Linking

Sd=LionLink(Sd,MeanBacLifed);
 S=LionLink(S,MeanBacLifeT);

% To do 2. close gaps and capture merging and splitting events. (Using Cost Matrix Gap Closing, merging, splitting.)


%% Combining + filtering spots

DeltaXcost=0.5;

Ilowerboundd=0;
IlowerboundT=0;

Sd=LionComBI(Sd,DeltaXcost,MeanBacLifed,Ilowerboundd);
% [Sd,dIntReduced]=LionCOMbee(Sd,DeltaXcost,MeanBacLifed,Ilowerboundd);

% dIntReduced is the COM method of spot combining, but accounting for one
% spot with that intensity compared to multiple spots with the same
% intensity.

S=LionComBI(S,DeltaXcost,MeanBacLifeT,IlowerboundT);


%% mean distance between spots of two channels (AFTER COMB and FILTERING)

[dd,ddweighted]=LionDistance(S,Sd,MeanBacLifeT,MeanBacLifed); 

%% Construct K for taking means of elements at same time point

SizeVectorT=[];
SizeVectord=[];

for i=1:Ncells
    
    SizeVectorT=[SizeVectorT size(S{i}.x,2)];
    SizeVectord=[SizeVectord size(Sd{i}.x,2)];
    
    NspotsT=size(S{i}.x,2);
    NspotsD=size(Sd{i}.x,2);
    
    for j=1:NspotsD
        
        Kdi{j}(:,i)=Sd{i}.x{j}(:,1);
        KdiFC{j}(:,i)=Sd{i}.x{j}(:,7); 
        
        Kdx{j}(:,i)=Sd{i}.x{j}(:,2);
        Kdy{j}(:,i)=Sd{i}.x{j}(:,4);      
        
    end
    for j=1:NspotsT
        Ki{j}(:,i)=S{i}.x{j}(:,1);
        KiFC{j}(:,i)=S{i}.x{j}(:,7);

        Kx{j}(:,i)=S{i}.x{j}(:,2);
        Ky{j}(:,i)=S{i}.x{j}(:,4);
    end
end


%% 
M=cell(Nspots,1); Mratio=zeros(Nspots,1); Mratiostd=zeros(Nspots,1);
Md=cell(Nspots,1); Mdratio=zeros(Nspots,1); Mdratiostd=zeros(Nspots,1);

for j=1:size(Ki,1)
        
        M{j}(:,1)=nanmean(Ki{j},2); % mean integrated intensity
        M{j}(:,2)=nanmean(Kx{j},2); % mean x position
        M{j}(:,3)=nanstd(Kx{j},2,1); % std x position
        M{j}(:,4)=nanmean(Ky{j},2); % mean y position
        M{j}(:,5)=nanstd(Ky{j},2,1); % std y position
        M{j}(:,6)=nanstd(Ki{j},2,1); % std intensities
        M{j}(:,7)=nanmean(KiFC{j},2); % mean full cell integrated intensity
        M{j}(:,8)=nanstd(KiFC{j},2,1); % std full cell integrated intensity
end
for j=1:size(Kdi,1)
        Md{j}(:,1)=nanmean(Kdi{j},2);
        Md{j}(:,2)=nanmean(Kdx{j},2);
        Md{j}(:,3)=nanstd(Kdx{j},2,1);
        Md{j}(:,4)=nanmean(Kdy{j},2);
        Md{j}(:,5)=nanstd(Kdy{j},2,1);
        Md{j}(:,6)=nanstd(Kdi{j},2,1);
        Md{j}(:,7)=nanmean(KdiFC{j},2);
        Md{j}(:,8)=nanstd(KdiFC{j},2,1);
end

%% Filtering

IupboundT=200000;
Iupboundd=150000;

n=Ncells+1;

T{n}.IntI=cell(1,max(SizeVectorT));
d{n}.IntI=cell(1,max(SizeVectord));

T{n}.FCI=cell(1,max(SizeVectorT));
d{n}.FCI=cell(1,max(SizeVectord));

T{n}.X=cell(1,max(SizeVectorT));
d{n}.X=cell(1,max(SizeVectord));

T{n}.Y=cell(1,max(SizeVectorT));
d{n}.Y=cell(1,max(SizeVectord));

for i=1:Ncells
    for j=1:size(S{i}.x,2) %that's Nspots for the T channel
    



% % Initiate Spot Integrated Intensity Cells
% 
% T{n}.IntI{j}=S{1}.x{j}(:,1);
% 
% % Initiate Full Cell Intensity Cells
% 
% T{n}.FCI{j}=S{1}.x{j}(:,7);
% 
% % Initiate normalised X Position Cells
% 
% T{n}.X{j}=S{1}.x{j}(:,2);
% 
% % Initiate normalised Y Position Cells
% 
% T{n}.Y{j}=S{1}.x{j}(:,4); 

    end
   for j=1:size(Sd{i}.x,2)
       % same for the d channel

           d{n}.IntI{j}=Sd{i}.x{j}(:,1);
           d{n}.FCI{j}=Sd{i}.x{j}(:,7);
           d{n}.X{j}=Sd{i}.x{j}(:,2);
           d{n}.Y{j}=Sd{i}.x{j}(:,4);
   end
   
end

% concatenate matrices to form 

for i=1:Ncells
    for j=1:size(S{i}.x,2)
        S{i}.x{j}(S{i}.x{j}>=IupboundT)=NaN;
        T{n}.FCI{j}=cat(1,T{n}.FCI{j},S{i}.x{j}(:,7));
        T{n}.X{j}=cat(1,T{n}.X{j},S{i}.x{j}(:,2));
        T{n}.Y{j}=cat(1,T{n}.Y{j},S{i}.x{j}(:,4));
        T{n}.IntI{j}=cat(1,T{n}.IntI{j},S{i}.x{j}(:,1));
    end
    for j=1:size(Sd{i}.x,2)
        Sd{i}.x{j}(Sd{i}.x{j}>=Iupboundd)=NaN;
        d{n}.FCI{j}=cat(1,d{n}.FCI{j},Sd{i}.x{j}(:,7));
        d{n}.X{j}=cat(1,d{n}.X{j},Sd{i}.x{j}(:,2));
        d{n}.Y{j}=cat(1,d{n}.Y{j},Sd{i}.x{j}(:,4));
        d{n}.IntI{j}=cat(1,d{n}.IntI{j},Sd{i}.x{j}(:,1));
    end
    
end

TotCellsStr=sprintf('Ncells = %d',Ncells);


for j=1:size(Ki,1)
MeanIntT{j}=mean(T{n}.IntI{j});
MeanIntFCT{j}=mean(T{n}.FCI{j});
StdIntT{j}=std(T{n}.IntI{j});
StdFCIntT{j}=std(T{n}.FCI{j});
MeanIntStrT{j}=sprintf('Mean = %g',MeanIntT{j});
StdIntStrT{j}=sprintf('std = %g',StdIntT{j});
end
for j=1:size(Kdi,1)
MeanIntd{j}=mean(d{n}.IntI{j});
MeanIntFCd{j}=mean(d{n}.FCI{j});
StdIntd{j}=std(d{n}.IntI{j});
StdFCIntd{j}=std(d{n}.FCI{j});
MeanIntStrd{j}=sprintf('Mean = %g',MeanIntd{j});
StdIntStrd{j}=sprintf('std = %g',StdIntd{j});
end

%% Reducing the variables for oversight.


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
