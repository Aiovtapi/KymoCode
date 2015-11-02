%% Load 2D Gaussian fit data and perform calculations.

clear all
close all
clc

% Tus - dif project
expno='001_DnaN_TUS_dif_30122014_TUSsignal';
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

BacLife=zeros(Ncells,1); BacLifed=zeros(Ncells,1);
MainPathTus=strcat(initval.basepath,'StacksLong/Tus/DataMULTI/');
MainPathdif=strcat(initval.basepath,'StacksLong/dif/DataMULTI/');


%% Load data
tic
for i=1:Ncells;
    T{i}=load(strcat(MainPathTus,num2str(i),'.mat'));
    d{i}=load(strcat(MainPathdif,num2str(i),'.mat'));
end

Nspots=T{i}.Nspots;
%% Calculations and Filtering
MeanIntT=cell(Nspots,1);MeanIntFCT=cell(Nspots,1);StdIntT=cell(Nspots,1);
StdFCIntT=cell(Nspots,1);MeanIntStrT=cell(Nspots,1);StdIntStrT=cell(Nspots,1);
StdIntd=cell(Nspots,1);StdFCIntd=cell(Nspots,1);MeanIntStrd=cell(Nspots,1);
StdIntStrd=cell(Nspots,1);

IupboundT=50000;
Iupboundd=500000;

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

T{n}.IntI{j}=T{1}.x{j}(:,6);
d{n}.IntI{j}=d{1}.x{j}(:,6);

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

T{n}.IntI{j}=cat(1,T{n}.IntI{j},T{i}.x{j}(:,6));
d{n}.IntI{j}=cat(1,d{n}.IntI{j},d{i}.x{j}(:,6));
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
    S{i}.x{j}(:,1)=imresize(T{i}.x{j}(:,6),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,7)=imresize(T{i}.x{j}(:,7),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,2)=imresize(T{i}.XNorm{j}(:,2),[MeanBacLifeT 1],'bilinear');
    S{i}.x{j}(:,4)=imresize(T{i}.XNorm{j}(:,4),[MeanBacLifeT 1],'bilinear');
    
    Sd{i}.x{j}(:,1)=imresize(d{i}.x{j}(:,6),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,7)=imresize(d{i}.x{j}(:,7),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,2)=imresize(d{i}.XNorm{j}(:,2),[MeanBacLifed 1],'bilinear');
    Sd{i}.x{j}(:,4)=imresize(d{i}.XNorm{j}(:,4),[MeanBacLifed 1],'bilinear');
    end
end

% Construct K for taking means of elements at same time point

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
        Kdy{j}(:,i)=Sd{i}.x(:,4);
        
    end
end

M=zeros(MeanBacLifeT,6);
Md=zeros(MeanBacLifeT,6);


for j=1:Nspots
        
        M{j}(:,1)=mean(Ki{j},2);
        M{j}(:,2)=mean(Kx{j},2);
        M{j}(:,3)=std(Kx{j},2);
        M{j}(:,4)=mean(Ky{j},2);
        M{j}(:,5)=std(Ky{j},2);
        M{j}(n,6)=std(Ki{j},2);
        M{j}(n,7)=mean(KiFC(n,:));
        M{j}(n,8)=std(KiFC(n,:));
        
end

for n=1:MeanBacLifed    
    Md(n,1)=mean(Kdi(n,:));
    MdR1(n,1)=mean(KdiR1(n,:));
    MdR2(n,1)=mean(KdiR2(n,:));
    MdR3(n,1)=mean(KdiR3(n,:));
    
    Md(n,2)=mean(Kdx(n,:));
    MdR1(n,2)=mean(KdxR1(n,:));
    MdR2(n,2)=mean(KdxR2(n,:));
    MdR3(n,2)=mean(KdxR3(n,:));
    
    Md(n,3)=std(Kdx(n,:));
    MdR1(n,3)=std(KdxR1(n,:));
    MdR2(n,3)=std(KdxR2(n,:));
    MdR3(n,3)=std(KdxR3(n,:));
    
    Md(n,4)=mean(Kdy(n,:));
    MdR1(n,4)=mean(KdyR1(n,:));
    MdR2(n,4)=mean(KdyR2(n,:));
    MdR3(n,4)=mean(KdyR3(n,:));
    
    Md(n,5)=std(Kdy(n,:));
    MdR1(n,5)=std(KdyR1(n,:));
    MdR2(n,5)=std(KdyR2(n,:));
    MdR3(n,5)=std(KdyR3(n,:));
    
    Md(n,6)=std(Kdi(n,:));
    MdR1(n,6)=std(KdiR1(n,:));
    MdR2(n,6)=std(KdiR2(n,:));
    MdR3(n,6)=std(KdiR3(n,:));
    
    Md(n,7)=mean(KdiFC(n,:));
    Md(n,8)=std(KdiFC(n,:));
end

Mratio=M(:,1)./M(:,7);
Mratiostd=std(Mratio);

Mdratio=Md(:,1)./Md(:,7);
Mdratiostd=std(Mdratio);

%% Estimation of Single Tus Protein 
SPItus=mean([(MeanIntT-MeanIntDR1) (MeanIntDR1-MeanIntDR2) ...
    (MeanIntDR2-MeanIntDR3)]);
%% Importing the real lifetime data from 'original' Bacpics
MainPathTusOri=strcat(initval.basepath,'Stacks/Tus/DataMULTI/');
MainPathdifOri=strcat(initval.basepath,'Stacks/dif/DataMULTI/');

for i=1:Ncells;
    Dori{i}=load(strcat(MainPathTusOri,num2str(i),'.mat'));
    dori{i}=load(strcat(MainPathdifOri,num2str(i),'.mat'));
    CelllifeT(i)=length(Dori{i}.XNorm(:,1));
    Celllifed(i)=length(dori{i}.XNorm(:,1));
end

MeanCellLifed=mean(Celllifed);
MeanCellLifeT=mean(CelllifeT);

toc
