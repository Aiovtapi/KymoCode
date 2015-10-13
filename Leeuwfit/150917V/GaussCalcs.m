%% Load 2D Gaussian fit data and perform calculations.

clear all
close all
clc

% Tus - dif project
expno='001_DnaN_TUS_dif_30122014_TUSsignal';
%initval=A001_Images_Set_Experiment(expno);      
initval.basepath='/Users/rleeuw/Work/Data/20141230_dnaN_dif_tus/dnaN_dif_tus_40msExpTime_5minAcqTimeYFP_30msExpCFP_002_C1/';

% oriZ - Dif project
%initval.basepath='/Users/rleeuw/Work/Data/OriZ-Dif_Results/';
%% Define variables

Ncells=18;

D=cell(Ncells+4,1);
d=cell(Ncells+4,1);
S=cell(Ncells,1);
Sd=cell(Ncells,1);
BacLife=zeros(Ncells,1); BacLifed=zeros(Ncells,1);
MainPathTus=strcat(initval.basepath,'StacksLong/Tus/DataMULTI/');
MainPathdif=strcat(initval.basepath,'StacksLong/dif/DataMULTI/');

%% Load data
tic
for i=1:Ncells;
    D{i}=load(strcat(MainPathTus,num2str(i),'.mat'));
    d{i}=load(strcat(MainPathdif,num2str(i),'.mat'));
end

%% Calculations and Filtering

IupboundD=50000;
Iupboundd=500000;

for i=1:Ncells
     % Filter out intensities by upper bound
     for j=1:length(D{i}.x(:,1))
    if D{i}.x(j,1)>=IupboundD
        D{i}.x(j,1)=0;
    end
    if D{i}.xR1(j,1)>=IupboundD
        D{i}.xR1(j,1)=0;
    end
    if D{i}.xR2(j,1)>=IupboundD
        D{i}.xR2(j,1)=0;
    end
    if D{i}.xR3(j,1)>=IupboundD
        D{i}.xR3(j,1)=0;
    end
    
    % Filter out integrated intensities by upper bound
    
    if D{i}.x(j,6)>=IupboundD
        D{i}.x(j,6)=0;
    end
    if D{i}.xR1(j,6)>=IupboundD
        D{i}.xR1(j,6)=0;
    end
    if D{i}.xR2(j,6)>=IupboundD
        D{i}.xR2(j,6)=0;
    end
    if D{i}.xR3(j,6)>=IupboundD
        D{i}.xR3(j,6)=0;
    end
     end
    % Filter out intensities by upper bound
     for j=1:length(d{i}.x(:,1))
     if d{i}.x(j,1) > Iupboundd
         d{i}.x(j,1)=0;
     end
     if d{i}.xR1(j,1) > Iupboundd
         d{i}.xR1(j,1)=0;
     end
     if d{i}.xR2(j,1) > Iupboundd
         d{i}.xR2(j,1)=0;
     end
     if d{i}.xR3(j,1) > Iupboundd
         d{i}.xR3(j,1)=0;
     end
     
     % Filter out integrated intensities by upper bound
     
     if d{i}.x(j,6)> Iupboundd
        d{i}.x(j,6)=0;
    end
    if d{i}.xR1(j,6)> Iupboundd
        d{i}.xR1(j,6)=0;
    end
    if d{i}.xR2(j,6)> Iupboundd
        d{i}.xR2(j,6)=0;
    end
    if d{i}.xR3(j,6)> Iupboundd
        d{i}.xR3(j,6)=0;
    end
     end
end


n=Ncells+1;

% Initiate Spot Intensity Cells

D{n}.I=D{1}.x(:,1);
D{n}.IR1=D{1}.xR1(:,1);
D{n}.IR2=D{1}.xR2(:,1);
D{n}.IR3=D{1}.xR3(:,1);

d{n}.I=d{1}.x(:,1);
d{n}.IR1=d{1}.xR1(:,1);
d{n}.IR2=d{1}.xR2(:,1);
d{n}.IR3=d{1}.xR3(:,1);

% Initiate Cell Intensity Cells

D{n}.FCI=D{1}.x(:,7);

d{n}.FCI=d{1}.x(:,7);

% Initiate X Position Cells

D{n}.X=D{1}.XNorm(:,2);
D{n}.XR1=D{1}.XNormR1(:,2);
D{n}.XR2=D{1}.XNormR2(:,2);
D{n}.XR3=D{1}.XNormR3(:,2);

d{n}.X=d{1}.XNorm(:,2);
d{n}.XR1=d{1}.XNormR1(:,2);
d{n}.XR2=d{1}.XNormR2(:,2);
d{n}.XR3=d{1}.XNormR3(:,2);

% Initiate Y Position Cells

D{n}.Y=D{1}.XNorm(:,4); 
D{n}.YR1=D{1}.XNormR1(:,4);
D{n}.YR2=D{1}.XNormR2(:,4);
D{n}.YR3=D{1}.XNormR3(:,4);

d{n}.Y=d{1}.XNorm(:,4);
d{n}.YR1=d{1}.XNormR1(:,4);
d{n}.YR2=d{1}.XNormR2(:,4);
d{n}.YR3=d{1}.XNormR3(:,4);

% Initiate Integrated Intensity Cells

D{n}.IntI=D{1}.x(:,6);
D{n}.IntIR1=D{1}.xR1(:,6);
D{n}.IntIR2=D{1}.xR2(:,6);
D{n}.IntIR3=D{1}.xR3(:,6);

d{n}.IntI=d{1}.x(:,6);
d{n}.IntIR1=d{1}.xR1(:,6);
d{n}.IntIR2=d{1}.xR2(:,6);
d{n}.IntIR3=d{1}.xR3(:,6);

for i=2:Ncells
D{n}.I=cat(1,D{n}.I,D{i}.x(:,1));
D{n}.IR1=cat(1,D{n}.IR1,D{i}.xR1(:,1));
D{n}.IR2=cat(1,D{n}.IR2,D{i}.xR2(:,1));
D{n}.IR3=cat(1,D{n}.IR3,D{i}.xR3(:,1));

D{n}.FCI=cat(1,D{n}.FCI,D{i}.x(:,7));

d{n}.FCI=cat(1,d{n}.FCI,d{i}.x(:,7));

d{n}.I=cat(1,d{n}.I,d{i}.x(:,1));
d{n}.IR1=cat(1,d{n}.IR1,d{i}.xR1(:,1));
d{n}.IR2=cat(1,d{n}.IR2,d{i}.xR2(:,1));
d{n}.IR3=cat(1,d{n}.IR3,d{i}.xR3(:,1));

D{n}.X=cat(1,D{n}.X,D{i}.XNorm(:,2));
D{n}.XR1=cat(1,D{n}.XR1,D{i}.XNormR1(:,2));
D{n}.XR2=cat(1,D{n}.XR2,D{i}.XNormR2(:,2));
D{n}.XR3=cat(1,D{n}.XR3,D{i}.XNormR3(:,2));

d{n}.X=cat(1,d{n}.X,d{i}.XNorm(:,2));
d{n}.XR1=cat(1,d{n}.XR1,d{i}.XNormR1(:,2));
d{n}.XR2=cat(1,d{n}.XR2,d{i}.XNormR2(:,2));
d{n}.XR3=cat(1,d{n}.XR3,d{i}.XNormR3(:,2));

D{n}.Y=cat(1,D{n}.Y,D{i}.XNorm(:,4));
D{n}.YR1=cat(1,D{n}.YR1,D{i}.XNormR1(:,4));
D{n}.YR2=cat(1,D{n}.YR2,D{i}.XNormR2(:,4));
D{n}.YR3=cat(1,D{n}.YR3,D{i}.XNormR3(:,4));

d{n}.Y=cat(1,d{n}.Y,d{i}.XNorm(:,4));
d{n}.YR1=cat(1,d{n}.YR1,d{i}.XNormR1(:,4));
d{n}.YR2=cat(1,d{n}.YR2,d{i}.XNormR2(:,4));
d{n}.YR3=cat(1,d{n}.YR3,d{i}.XNormR3(:,4));

D{n}.IntI=cat(1,D{n}.IntI,D{i}.x(:,6));
D{n}.IntIR1=cat(1,D{n}.IntIR1,D{i}.xR1(:,6));
D{n}.IntIR2=cat(1,D{n}.IntIR2,D{i}.xR2(:,6));
D{n}.IntIR3=cat(1,D{n}.IntIR3,D{i}.xR3(:,6));

d{n}.IntI=cat(1,d{n}.IntI,d{i}.x(:,6));
d{n}.IntIR1=cat(1,d{n}.IntIR1,d{i}.xR1(:,6));
d{n}.IntIR2=cat(1,d{n}.IntIR2,d{i}.xR2(:,6));
d{n}.IntIR3=cat(1,d{n}.IntIR3,d{i}.xR3(:,6));
end

%% Spot tracking calculations -- switch positions between spots given criterium
% exchange the spot positions with ones closest to previous valuesSp

for i=1:Ncells
    for j=2:length(D{i}.x(:,1));
        
        CritValue=0.2;
        
Diffvectorx=[abs(D{i}.XNorm(j,2)-D{i}.XNorm(j-1,2)) ...
            abs(D{i}.XNormR1(j,2)-D{i}.XNorm(j-1,2)) ...
            abs(D{i}.XNormR2(j,2)-D{i}.XNorm(j-1,2)) ...
            abs(D{i}.XNormR3(j,2)-D{i}.XNorm(j-1,2))];
        
[Vdiff,idx]=min(Diffvectorx);

if Vdiff<=CritValue
    DUM=D{i}.XNorm(j,2);
    if idx==2;
    D{i}.XNorm(j,2)=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=DUM;
    elseif idx==3;
    D{i}.XNorm(j,2)=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=DUM;
    elseif idx==4;
    D{i}.XNorm(j,2)=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=DUM;
    end
end

DiffvectorxR1=[abs(D{i}.XNorm(j,2)-D{i}.XNormR1(j-1,2)) ...
            abs(D{i}.XNormR1(j,2)-D{i}.XNormR1(j-1,2)) ...
            abs(D{i}.XNormR2(j,2)-D{i}.XNormR1(j-1,2)) ...
            abs(D{i}.XNormR3(j,2)-D{i}.XNormR1(j-1,2))];
        
[VdiffR1,idxR1]=min(DiffvectorxR1);


if VdiffR1<=CritValue
    DUM=D{i}.XNormR1(j,2);
    if idxR1==1;
    D{i}.XNormR1(j,2)=D{i}.XNorm(j,2);
    D{i}.XNorm(j,2)=DUM;
    elseif idxR1==3;
    D{i}.XNormR1(j,2)=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=DUM;
    elseif idxR1==4;
    D{i}.XNormR1(j,2)=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=DUM;
    end
end
  
DiffvectorxR2=[abs(D{i}.XNorm(j,2)-D{i}.XNormR2(j-1,2)) ...
            abs(D{i}.XNormR1(j,2)-D{i}.XNormR2(j-1,2)) ...
            abs(D{i}.XNormR2(j,2)-D{i}.XNormR2(j-1,2)) ...
            abs(D{i}.XNormR3(j,2)-D{i}.XNormR2(j-1,2))];
        
[VdiffR2,idxR2]=min(DiffvectorxR2);

if VdiffR2<=CritValue
    DUM=D{i}.XNormR2(j,2);
    if idxR2==1;
    D{i}.XNormR2(j,2)=D{i}.XNorm(j,2);
    D{i}.XNorm(j,2)=DUM;
    elseif idxR2==2;
    D{i}.XNormR2(j,2)=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=DUM;
    elseif idxR2==4;
    D{i}.XNormR2(j,2)=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=DUM;
    end
end

DiffvectorxR3=[abs(D{i}.XNorm(j,2)-D{i}.XNormR3(j-1,2)) ...
            abs(D{i}.XNormR1(j,2)-D{i}.XNormR3(j-1,2)) ...
            abs(D{i}.XNormR2(j,2)-D{i}.XNormR3(j-1,2)) ...
            abs(D{i}.XNormR3(j,2)-D{i}.XNormR3(j-1,2))];
        
[VdiffR3,idxR3]=min(DiffvectorxR3);

if VdiffR3<=CritValue
    DUM=D{i}.XNormR3(j,2);
    if idxR3==1;
    D{i}.XNormR3(j,2)=D{i}.XNorm(j,2);
    D{i}.XNorm(j,2)=DUM;
    elseif idxR3==2;
    D{i}.XNormR3(j,2)=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=DUM;
    elseif idxR3==3;
    D{i}.XNormR3(j,2)=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=DUM;
    end
end

    end
end

%% Some Calculations
TotCellsStr=sprintf('Ncells = %d',Ncells);

MeanIntD=mean(D{n}.IntI);
MeanIntDR1=mean(D{n}.IntIR1);
MeanIntDR2=mean(D{n}.IntIR2);
MeanIntDR3=mean(D{n}.IntIR3);

MeanIntFCD=mean(D{n}.FCI);

StdIntD=std(D{n}.IntI);
StdIntDR1=std(D{n}.IntIR1);
StdIntDR2=std(D{n}.IntIR2);
StdIntDR3=std(D{n}.IntIR3);

StdFCIntD=std(D{n}.FCI);

MeanIntStrD=sprintf('Mean = %g',MeanIntD);
MeanIntStrDR1=sprintf('Mean = %g',MeanIntDR1);
MeanIntStrDR2=sprintf('Mean = %g',MeanIntDR2);
MeanIntStrDR3=sprintf('Mean = %g',MeanIntDR3);

StdIntStrD=sprintf('std = %g',StdIntD);
StdIntStrDR1=sprintf('std = %g',StdIntDR1);
StdIntStrDR2=sprintf('std = %g',StdIntDR2);
StdIntStrDR3=sprintf('std = %g',StdIntDR3);

MeanIntd=mean(d{n}.IntI);
MeanIntdR1=mean(d{n}.IntIR1);
MeanIntdR2=mean(d{n}.IntIR2);
MeanIntdR3=mean(d{n}.IntIR3);

MeanIntFCd=mean(d{n}.FCI);

StdIntd=std(d{n}.IntI);
StdIntdR1=std(d{n}.IntIR1);
StdIntdR2=std(d{n}.IntIR2);
StdIntdR3=std(d{n}.IntIR3);

StdFCIntd=std(d{n}.FCI);

MeanIntStrd=sprintf('Mean = %g',MeanIntd);
MeanIntStrdR1=sprintf('Mean = %g',MeanIntdR1);
MeanIntStrdR2=sprintf('Mean = %g',MeanIntdR2);
MeanIntStrdR3=sprintf('Mean = %g',MeanIntdR3);

StdIntStrd=sprintf('std = %g',StdIntd);
StdIntStrdR1=sprintf('std = %g',StdIntdR1);
StdIntStrdR2=sprintf('std = %g',StdIntdR2);
StdIntStrdR3=sprintf('std = %g',StdIntdR3);

for i=1:Ncells
    BacLife(i)=length(D{i}.x(:,1));
    BacLifed(i)=length(d{i}.x(:,1));
end

MaxBacLife=30;
MaxBacLifed=max(BacLifed);

Ki=zeros(MaxBacLife,Ncells);
KiFC=zeros(MaxBacLife,Ncells);
KiR1=zeros(MaxBacLife,Ncells);
KiR2=zeros(MaxBacLife,Ncells);
KiR3=zeros(MaxBacLife,Ncells);

Kx=zeros(MaxBacLife,Ncells);
KxR1=zeros(MaxBacLife,Ncells);
KxR2=zeros(MaxBacLife,Ncells);
KxR3=zeros(MaxBacLife,Ncells);

Ky=zeros(MaxBacLife,Ncells);
KyR1=zeros(MaxBacLife,Ncells);
KyR2=zeros(MaxBacLife,Ncells);
KyR3=zeros(MaxBacLife,Ncells);

Kdi=zeros(MaxBacLife,Ncells);
KdiFC=zeros(MaxBacLife,Ncells);
KdiR1=zeros(MaxBacLifed,Ncells);
KdiR2=zeros(MaxBacLifed,Ncells);
KdiR3=zeros(MaxBacLifed,Ncells);

Kdx=zeros(MaxBacLife,Ncells);
KdxR1=zeros(MaxBacLifed,Ncells);
KdxR2=zeros(MaxBacLifed,Ncells);
KdxR3=zeros(MaxBacLifed,Ncells);

Kdy=zeros(MaxBacLife,Ncells);
KdyR1=zeros(MaxBacLifed,Ncells);
KdyR2=zeros(MaxBacLifed,Ncells);
KdyR3=zeros(MaxBacLifed,Ncells);

for i=1:Ncells
    S{i}.x(:,1)=imresize(D{i}.x(:,6),[MaxBacLife 1],'bilinear');
    S{i}.xR1(:,1)=imresize(D{i}.xR1(:,6),[MaxBacLife 1],'bilinear');
    S{i}.xR2(:,1)=imresize(D{i}.xR2(:,6),[MaxBacLife 1],'bilinear');
    S{i}.xR3(:,1)=imresize(D{i}.xR3(:,6),[MaxBacLife 1],'bilinear');
    
    S{i}.x(:,7)=imresize(D{i}.x(:,7),[MaxBacLife 1],'bilinear');
    
    S{i}.x(:,2)=imresize(D{i}.XNorm(:,2),[MaxBacLife 1],'bilinear');
    S{i}.xR1(:,2)=imresize(D{i}.XNormR1(:,2),[MaxBacLife 1],'bilinear');
    S{i}.xR2(:,2)=imresize(D{i}.XNormR2(:,2),[MaxBacLife 1],'bilinear');
    S{i}.xR3(:,2)=imresize(D{i}.XNormR3(:,2),[MaxBacLife 1],'bilinear');
    
    S{i}.x(:,4)=imresize(D{i}.XNorm(:,4),[MaxBacLife 1],'bilinear');
    S{i}.xR1(:,4)=imresize(D{i}.XNormR1(:,4),[MaxBacLife 1],'bilinear');
    S{i}.xR2(:,4)=imresize(D{i}.XNormR2(:,4),[MaxBacLife 1],'bilinear');
    S{i}.xR3(:,4)=imresize(D{i}.XNormR3(:,4),[MaxBacLife 1],'bilinear');
    
    Sd{i}.x(:,1)=imresize(d{i}.x(:,6),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR1(:,1)=imresize(d{i}.xR1(:,6),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR2(:,1)=imresize(d{i}.xR2(:,6),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR3(:,1)=imresize(d{i}.xR3(:,6),[MaxBacLifed 1],'bilinear');
    
    Sd{i}.x(:,7)=imresize(d{i}.x(:,7),[MaxBacLifed 1],'bilinear');
    
    Sd{i}.x(:,2)=imresize(d{i}.XNorm(:,2),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR1(:,2)=imresize(d{i}.XNormR1(:,2),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR2(:,2)=imresize(d{i}.XNormR2(:,2),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR3(:,2)=imresize(d{i}.XNormR3(:,2),[MaxBacLifed 1],'bilinear');
    
    Sd{i}.x(:,4)=imresize(d{i}.XNorm(:,4),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR1(:,4)=imresize(d{i}.XNormR1(:,4),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR2(:,4)=imresize(d{i}.XNormR2(:,4),[MaxBacLifed 1],'bilinear');
    Sd{i}.xR3(:,4)=imresize(d{i}.XNormR3(:,4),[MaxBacLifed 1],'bilinear');
    
end

%Construct K for taking means of elements at same time point
for j=1:MaxBacLife
    for i=1:Ncells
        Ki(j,i)=S{i}.x(j,1);
        KiFC(j,i)=S{i}.x(j,7);
        KiR1(j,i)=S{i}.xR1(j,1);
        KiR2(j,i)=S{i}.xR2(j,1);
        KiR3(j,i)=S{i}.xR3(j,1);

        Kx(j,i)=S{i}.x(j,2);
        KxR1(j,i)=S{i}.xR1(j,2);
        KxR2(j,i)=S{i}.xR2(j,2);
        KxR3(j,i)=S{i}.xR3(j,2);
        
        Ky(j,i)=S{i}.x(j,4);
        KyR1(j,i)=S{i}.xR1(j,4);
        KyR2(j,i)=S{i}.xR2(j,4);
        KyR3(j,i)=S{i}.xR3(j,4);
    end
end

for j=1:MaxBacLifed
    for i=1:Ncells
        Kdi(j,i)=Sd{i}.x(j,1);
        KdiFC(j,i)=Sd{i}.x(j,7);
        KdiR1(j,i)=Sd{i}.xR1(j,1);
        KdiR2(j,i)=Sd{i}.xR2(j,1);
        KdiR3(j,i)=Sd{i}.xR3(j,1);
        
        Kdx(j,i)=Sd{i}.x(j,2);
        KdxR1(j,i)=Sd{i}.xR1(j,2);
        KdxR2(j,i)=Sd{i}.xR2(j,2);
        KdxR3(j,i)=Sd{i}.xR3(j,2);
        
        Kdy(j,i)=Sd{i}.x(j,4);
        KdyR1(j,i)=Sd{i}.xR1(j,4);
        KdyR2(j,i)=Sd{i}.xR2(j,4);
        KdyR3(j,i)=Sd{i}.xR3(j,4);
        
    end
end

M=zeros(MaxBacLife,6);
MR1=zeros(MaxBacLife,6);
MR2=zeros(MaxBacLife,6);
MR3=zeros(MaxBacLife,6);

Md=zeros(MaxBacLife,6);
MdR1=zeros(MaxBacLife,6);
MdR2=zeros(MaxBacLife,6);
MdR3=zeros(MaxBacLife,6);

for j=1:MaxBacLife
    M(j,1)=mean(Ki(j,:));
    MR1(j,1)=mean(KiR1(j,:));
    MR2(j,1)=mean(KiR2(j,:));
    MR3(j,1)=mean(KiR3(j,:));
    
    M(j,2)=mean(Kx(j,:));
    MR1(j,2)=mean(KxR1(j,:));
    MR2(j,2)=mean(KxR2(j,:));
    MR3(j,2)=mean(KxR3(j,:));
    
    M(j,3)=std(Kx(j,:));
    MR1(j,3)=std(KxR1(j,:));
    MR2(j,3)=std(KxR2(j,:));
    MR3(j,3)=std(KxR3(j,:));
    
    M(j,4)=mean(Ky(j,:));
    MR1(j,4)=mean(KyR1(j,:));
    MR2(j,4)=mean(KyR2(j,:));
    MR3(j,4)=mean(KyR3(j,:));
    
    M(j,5)=std(Ky(j,:));
    MR1(j,5)=std(KyR1(j,:));
    MR2(j,5)=std(KyR2(j,:));
    MR3(j,5)=std(KyR3(j,:));
    
    M(j,6)=std(Ki(j,:));
    MR1(j,6)=std(KiR1(j,:));
    MR2(j,6)=std(KiR2(j,:));
    MR3(j,6)=std(KiR3(j,:));
    
    M(j,7)=mean(KiFC(j,:));
    M(j,8)=std(KiFC(j,:));
end

for j=1:MaxBacLifed    
    Md(j,1)=mean(Kdi(j,:));
    MdR1(j,1)=mean(KdiR1(j,:));
    MdR2(j,1)=mean(KdiR2(j,:));
    MdR3(j,1)=mean(KdiR3(j,:));
    
    Md(j,2)=mean(Kdx(j,:));
    MdR1(j,2)=mean(KdxR1(j,:));
    MdR2(j,2)=mean(KdxR2(j,:));
    MdR3(j,2)=mean(KdxR3(j,:));
    
    Md(j,3)=std(Kdx(j,:));
    MdR1(j,3)=std(KdxR1(j,:));
    MdR2(j,3)=std(KdxR2(j,:));
    MdR3(j,3)=std(KdxR3(j,:));
    
    Md(j,4)=mean(Kdy(j,:));
    MdR1(j,4)=mean(KdyR1(j,:));
    MdR2(j,4)=mean(KdyR2(j,:));
    MdR3(j,4)=mean(KdyR3(j,:));
    
    Md(j,5)=std(Kdy(j,:));
    MdR1(j,5)=std(KdyR1(j,:));
    MdR2(j,5)=std(KdyR2(j,:));
    MdR3(j,5)=std(KdyR3(j,:));
    
    Md(j,6)=std(Kdi(j,:));
    MdR1(j,6)=std(KdiR1(j,:));
    MdR2(j,6)=std(KdiR2(j,:));
    MdR3(j,6)=std(KdiR3(j,:));
    
    Md(j,7)=mean(KdiFC(j,:));
    Md(j,8)=std(KdiFC(j,:));
end

Mratio=M(:,1)./M(:,7);
Mratiostd=std(Mratio);

Mdratio=Md(:,1)./Md(:,7);
Mdratiostd=std(Mdratio);

%% Estimation of Single Tus Protein 
SPItus=mean([(MeanIntD-MeanIntDR1) (MeanIntDR1-MeanIntDR2) ...
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
