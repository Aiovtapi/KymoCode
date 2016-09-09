%% Load 2D Gaussian fit data and perform calculations.

clear all
close all
clc
 
% oriZ - Dif project
%initval.basepath='/Users/rleeuw/Work/Data/OriZ-Dif_Results/';
%% Define variables

Ncells=31;

T=cell(Ncells+1,1);
d=cell(Ncells+1,1);
S=cell(Ncells,1);
Sd=cell(Ncells,1);
SnN=cell(Ncells,1);
SdnN=cell(Ncells,1);
framesd=zeros(Ncells,1);
framesT=zeros(Ncells,1);

% Tus - dif project

%  expno='001_DnaN_TUS_dif_30122014_difsignal';
%  initval=A001_Images_Set_Experiment(expno); %define your paths and files


% Tus Manual Labor!
BacLife=zeros(Ncells,1); BacLifed=zeros(Ncells,1);

lionval.channr=7;
lionval.viewchan='CFP';
lionval.viewbac=1:22;

exp='Roy_MM_Tus_dif';

[lionval]=LionDefine(exp,lionval);

%% Load data
tic

%SIMULATION FOR MSDs

for i=1:Ncells;
    T{i}=load(strcat(lionval.MainPathTus,'cell',num2str(i,'%03.0f'),'.mat'));
    framesT(i)=length(T{i}.x{1}(:,1));
    d{i}=load(strcat(lionval.MainPathdif,'cell',num2str(i,'%03.0f'),'.mat'));
    framesd(i)=length(d{i}.x{1}(:,1));
end

NspotsT=T{i}.NSpots;
NspotsD=d{i}.NSpots;
Npointsd=sum(framesd); 
NpointsT=sum(framesT); 
%% Calculations and Filtering

MeanIntT=cell(NspotsT,1);MeanIntFCT=cell(NspotsT,1);StdIntT=cell(NspotsT,1);
StdFCIntT=cell(NspotsT,1);MeanIntStrT=cell(NspotsT,1);StdIntStrT=cell(NspotsT,1);
StdIntd=cell(NspotsD,1);StdFCIntd=cell(NspotsD,1);MeanIntStrd=cell(NspotsD,1);
StdIntStrd=cell(NspotsD,1);

n=Ncells+1;


%% Some Calculations


for i=1:Ncells
    BacLife(i)=length(T{i}.x{1}(:,1));
    BacLifed(i)=length(d{i}.x{1}(:,1));
end

Ki=cell(NspotsT,1);
KiFC=cell(NspotsT,1);
Kx=cell(NspotsT,1);
Ky=cell(NspotsT,1);

Kdi=cell(NspotsD,1);
KdiFC=cell(NspotsD,1);
Kdx=cell(NspotsD,1);
Kdy=cell(NspotsD,1);

MeanBacLifeT=round(mean(BacLife));
MeanBacLifed=round(mean(BacLifed));

% I can't resize images now?

for i=1:Ncells
    for j=1:size(d{i}.x,2) 
    SdS{i}.x{j}(:,1)=imresize(d{i}.x{j}(:,8),[MeanBacLifed 1],'bilinear');
    SdS{i}.x{j}(:,7)=imresize(d{i}.x{j}(:,7),[MeanBacLifed 1],'bilinear');
    SdS{i}.x{j}(:,3)=imresize(d{i}.x{j}(:,3),[MeanBacLifed 1],'bilinear');
    SdS{i}.x{j}(:,5)=imresize(d{i}.x{j}(:,5),[MeanBacLifed 1],'bilinear');
    SdS{i}.x{j}(:,2)=imresize(d{i}.XNorm{j}(:,2),[MeanBacLifed 1],'bilinear');
    SdS{i}.x{j}(:,4)=imresize(d{i}.XNorm{j}(:,4),[MeanBacLifed 1],'bilinear');
    SdS{i}.x{j}(:,8)=imresize(d{i}.pixels{j}(1,:)',[MeanBacLifed 1],'bilinear');
    end
    
    for j=1:size(T{i}.x,2)
    SS{i}.x{j}(:,1)=imresize(T{i}.x{j}(:,8),[MeanBacLifeT 1],'bilinear');
    SS{i}.x{j}(:,7)=imresize(T{i}.x{j}(:,7),[MeanBacLifeT 1],'bilinear');
    SS{i}.x{j}(:,3)=imresize(T{i}.x{j}(:,3),[MeanBacLifeT 1],'bilinear');
    SS{i}.x{j}(:,3)=imresize(T{i}.x{j}(:,5),[MeanBacLifeT 1],'bilinear');
    SS{i}.x{j}(:,2)=imresize(T{i}.XNorm{j}(:,2),[MeanBacLifeT 1],'bilinear');
    SS{i}.x{j}(:,4)=imresize(T{i}.XNorm{j}(:,4),[MeanBacLifeT 1],'bilinear');
    SS{i}.x{j}(:,8)=imresize(T{i}.pixels{j}(1,:)',[MeanBacLifeT 1],'bilinear');
    end
end

for i=1:Ncells
    for j=1:size(d{i}.x,2)
      Sd{i}.x{j}(:,1)=d{i}.x{j}(:,8);
      Sd{i}.x{j}(:,2)=d{i}.XNorm{j}(:,2);
      Sd{i}.x{j}(:,3)=d{i}.x{j}(:,3);
      Sd{i}.x{j}(:,4)=d{i}.XNorm{j}(:,4);
      Sd{i}.x{j}(:,5)=d{i}.x{j}(:,5);
      Sd{i}.x{j}(:,7)=d{i}.x{j}(:,7);
      
      SdnN{i}.x{j}(:,1)=d{i}.x{j}(:,8);
      SdnN{i}.x{j}(:,2)=d{i}.x{j}(:,2);
      SdnN{i}.x{j}(:,3)=d{i}.x{j}(:,3);
      SdnN{i}.x{j}(:,4)=d{i}.x{j}(:,4);
      SdnN{i}.x{j}(:,5)=d{i}.x{j}(:,5);
      SdnN{i}.x{j}(:,7)=d{i}.x{j}(:,7);
    end
    
    for j=1:size(T{i}.x,2)
      S{i}.x{j}(:,1)=T{i}.x{j}(:,8);
      S{i}.x{j}(:,2)=T{i}.XNorm{j}(:,2);
      S{i}.x{j}(:,3)=T{i}.x{j}(:,3);
      S{i}.x{j}(:,4)=T{i}.XNorm{j}(:,4);
      S{i}.x{j}(:,5)=T{i}.x{j}(:,5);
      S{i}.x{j}(:,7)=T{i}.x{j}(:,7);
      
      SnN{i}.x{j}(:,1)=T{i}.x{j}(:,8);
      SnN{i}.x{j}(:,2)=T{i}.x{j}(:,2);
      SnN{i}.x{j}(:,3)=T{i}.x{j}(:,3);
      SnN{i}.x{j}(:,4)=T{i}.x{j}(:,4);
      SnN{i}.x{j}(:,5)=T{i}.x{j}(:,5);
      SnN{i}.x{j}(:,7)=T{i}.x{j}(:,7);
    end
end


%% Combining + intensity filtering spots

Ilowerboundd=1000;
IlowerboundT=500;

pixel_separation_constant=4;

SdS=LionComBI(SdS,d,pixel_separation_constant,Ilowerboundd);
SS=LionComBI(SS,T,pixel_separation_constant,IlowerboundT);

% % [Sd,dIntReduced]=LionCOMbee(Sd,DeltaXcost,MeanBacLifed,Ilowerboundd);
% 
% % dIntReduced is the COM method of spot combining, but accounting for one
% % spot with that intensity compared to multiple spots with the same
% % intensity.
% 
% S=LionComBI(S,DeltaXcost,MeanBacLifeT,IlowerboundT);
% SS=LionComBIS(SS,DeltaXcost,MeanBacLifeT,IlowerboundT);

%% Spot tracking + linking algorithms
%1. Linking

S=LionLink(S);
%  SnN=LionLink(SnN);

% To do 2. close gaps and capture merging and splitting events. (Using Cost Matrix Gap Closing, merging, splitting.)

%% mean distance between spots of two channels (AFTER COMB and FILTERING)

[dd,ddweighted]=LionDistance(S,Sd,MeanBacLifeT,MeanBacLifed); 

%% Construct K for taking means of elements at same time point

SizeVectorT=[];
SizeVectord=[];

for i=1:Ncells
    
    SizeVectorT=[SizeVectorT size(S{i}.x,2)];
    SizeVectord=[SizeVectord size(Sd{i}.x,2)];
    
    NspotsT=size(S{i}.x,2);
    NspotsTS=size(SS{i}.x,2);

    NspotsD=size(Sd{i}.x,2);
    NspotsDS=size(SdS{i}.x,2);
    
    for j=1:NspotsDS
        
        Kdi{j}(:,i)=SdS{i}.x{j}(:,1);
        KdiFC{j}(:,i)=SdS{i}.x{j}(:,7); 
        
        Kdx{j}(:,i)=SdS{i}.x{j}(:,2);
        Kdy{j}(:,i)=SdS{i}.x{j}(:,4);      
        
        Kdpx{j}(:,i)=SdS{i}.x{j}(:,8);
    end
    for j=1:NspotsTS
        Ki{j}(:,i)=SS{i}.x{j}(:,1);
        KiFC{j}(:,i)=SS{i}.x{j}(:,7);

        Kx{j}(:,i)=SS{i}.x{j}(:,2);
        Ky{j}(:,i)=SS{i}.x{j}(:,4);
        
        Kpx{j}(:,i)=SS{i}.x{j}(:,8);
    end
end


%% 
M=cell(NspotsT,1); Mratio=zeros(NspotsT,1); Mratiostd=zeros(NspotsT,1);
Md=cell(NspotsD,1); Mdratio=zeros(NspotsD,1); Mdratiostd=zeros(NspotsD,1);

% for j=1:NspotsT
%         
%         M{j}(:,1)=nanmean(Ki{j},2); % mean integrated intensity
%         M{j}(:,2)=nanmean(Kx{j},2); % mean x position
%         M{j}(:,3)=nanstd(Kx{j},1); % std x position
%         M{j}(:,4)=nanmean(Ky{j},2); % mean y position
%         M{j}(:,5)=nanstd(Ky{j},1); % std y position
%         M{j}(:,6)=nanstd(Ki{j},1); % std intensities
%         M{j}(:,7)=nanmean(KiFC{j},2); % mean full cell integrated intensity
%         M{j}(:,8)=nanstd(KiFC{j},1); % std full cell integrated intensity
% end
% for j=1:NspotsD
%         Md{j}(:,1)=nanmean(Kdi{j},2);
%         Md{j}(:,2)=nanmean(Kdx{j},2);
%         Md{j}(:,3)=nanstd(Kdx{j},1);
%         Md{j}(:,4)=nanmean(Kdy{j},2);
%         Md{j}(:,5)=nanstd(Kdy{j},1);
%         Md{j}(:,6)=nanstd(Kdi{j},1);
%         Md{j}(:,7)=nanmean(KdiFC{j},2);
%         Md{j}(:,8)=nanstd(KdiFC{j},1);
% end

%"I" is the integrated intensity of the spot.
[dif.I, dif.Istd]=LionRowMS(Kdi);
[dif.x, dif.xstd]=LionRowMS(Kdx);
[dif.y, dif.ystd]=LionRowMS(Kdy);
[dif.FC, dif.FCstd]=LionRowMS(KdiFC);

[Tus.I, Tus.Istd]=LionRowMS(Ki);
[Tus.x, Tus.xstd]=LionRowMS(Kx);
[Tus.y, Tus.ystd]=LionRowMS(Ky);
[Tus.FC, Tus.FCstd]=LionRowMS(KiFC);

[Tus.px, Tus.pxstd]=LionRowMS(Kpx');
[dif.px, dif.pxstd]=LionRowMS(Kdpx');

%% Filtering

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
        SpotNumbersd(i)=size(Sd{i}.x,2);
        SpotMaxNumd=max(SpotNumbersd);
        
        SpotNumbersT(i)=size(S{i}.x,2);
        SpotMaxNumT=max(SpotNumbersT);
end


for i=1:Ncells
    for j=1:size(S{i}.x,2) 
        %that's Nspots for the T channel
        %(But this does not line up well yet!)
           T{n}.IntI{j}=S{i}.x{j}(:,1);
           T{n}.FCI{j}=S{i}.x{j}(:,7);
           T{n}.X{j}=S{i}.x{j}(:,2);
           T{n}.Y{j}=S{i}.x{j}(:,4);

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

XdumVecd=cell(Ncells,1);

for i=1:Ncells
    for j=1:size(S{i}.x,2)
%         S{i}.x{j}(S{i}.x{j}>=IupboundT)=NaN;
        T{n}.FCI{j}=cat(1,T{n}.FCI{j},S{i}.x{j}(:,7));
        T{n}.X{j}=cat(1,T{n}.X{j},S{i}.x{j}(:,2));
        T{n}.Y{j}=cat(1,T{n}.Y{j},S{i}.x{j}(:,4));
        T{n}.IntI{j}=cat(1,T{n}.IntI{j},S{i}.x{j}(:,1));
    end
    for j=1:size(Sd{i}.x,2)
%         Sd{i}.x{j}(Sd{i}.x{j}>=Iupboundd)=NaN;
        d{n}.FCI{j}=cat(1,d{n}.FCI{j},Sd{i}.x{j}(:,7));
        d{n}.X{j}=cat(1,d{n}.X{j},Sd{i}.x{j}(:,2));
        d{n}.Y{j}=cat(1,d{n}.Y{j},Sd{i}.x{j}(:,4));
        d{n}.IntI{j}=cat(1,d{n}.IntI{j},Sd{i}.x{j}(:,1));
        
        
        XdumVecd{i}(:,j)=Sd{i}.x{j}(:,2);
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

dXmat=cell(Ncells,1);
dYmat=cell(Ncells,1);
dImat=cell(Ncells,1);

TXmat=cell(Ncells,1);
TYmat=cell(Ncells,1);
TImat=cell(Ncells,1);


for i=1:Ncells
    for j=1:size(S{i}.x,2)
        
            TXmat{i}(:,j)=S{i}.x{j}(:,2);
            TImat{i}(:,j)=S{i}.x{j}(:,1);
            TYmat{i}(:,j)=S{i}.x{j}(:,4);
    end
    
    for j=1:size(Sd{i}.x,2)
                    
            dXmat{i}(:,j)=Sd{i}.x{j}(:,2);
            dImat{i}(:,j)=Sd{i}.x{j}(:,1);
            dYmat{i}(:,j)=Sd{i}.x{j}(:,4);
            
    end
    
    dXmat{i}(isnan(dXmat{i}))=0;
    
    for k=1:length(dXmat{i}(:,1))
    IdXmat{k,i}=find(dXmat{i}(k,:));
    end
    for k=1:length(TXmat{i}(:,1))
    ITXmat{k,i}=find(TXmat{i}(k,:));
    end
    
end


%% additional things
Tus.totalspotintensity=zeros(length(Tus.x{1}),1);
dif.totalspotintensity=zeros(length(dif.x{1}),1);

Tus.totalspotintensitystd=zeros(length(Tus.x{1}),1);
dif.totalspotintensitystd=zeros(length(dif.x{1}),1);

Tus.totalcellintensity=zeros(length(Tus.x{1}),1);
dif.totalcellintensity=zeros(length(dif.x{1}),1);

Tus.totalcellintensitystd=zeros(length(Tus.x{1}),1);
dif.totalcellintensitystd=zeros(length(dif.x{1}),1);

Tus.activespots=zeros(length(Tus.I{1}),1);
dif.activespots=zeros(length(dif.I{1}),1);

for i=1:size(Tus.x,2) %spotnum
Tus.totalspotintensity=Tus.totalspotintensity+Tus.I{i};

Tus.totalspotintensitystd=Tus.totalspotintensitystd+Tus.Istd{i}.^2;

Tus.totalcellintensity=Tus.totalcellintensity+Tus.FC{i};

    for j=1:size(Tus.I{i},1)
        F=find(Tus.I{i}(j,1)>IlowerboundT);
        if isempty(F)
            F=0;
        end
    Tus.activespots(j)=Tus.activespots(j)+F;
    end
end

Tus.totalspotintensitystd=sqrt(Tus.totalspotintensity);

for m=1:size(dif.x,2)
dif.totalspotintensity=dif.totalspotintensity+dif.I{m};
dif.totalspotintensitystd=dif.totalspotintensitystd+dif.Istd{m}.^2;

dif.totalcellintensity=dif.totalcellintensity+dif.FC{m};
    for j=1:size(dif.I{m},1)
        F=find(dif.I{m}(j,1)>Ilowerboundd);
        if isempty(F)
            F=0;
        end
    dif.activespots(j)=dif.activespots(j)+F;
    end
end

dif.totalspotintensitystd=sqrt(dif.totalspotintensity);


%% Importing the real lifetime data from 'original' Bacpics
% MainPathTusOri=strcat(initval.basepath,'Stacks/Tus/DataMULTI/');
% MainPathdifOri=strcat(initval.basepath,'Stacks/dif/DataMULTI/');
% 
% for i=1:Ncells;
%     Tori{i}=load(strcat(MainPathTusOri,num2str(i),'.mat'));
%     dori{i}=load(strcat(MainPathdifOri,num2str(i),'.mat'));
%     CelllifeT(i)=length(Tori{i}.XNorm(:,1));
%     Celllifed(i)=length(dori{i}.XNorm(:,1));
% end
% 
% MeanCellLifed=mean(Celllifed);
% MeanCellLifeT=mean(CelllifeT);

%% Cell Length + mean growth
for i=1:Ncells
    for t=1:size(d{i}.ydatacrpdR1,1)
        dif.length{i}(t,1)=size(d{i}.ydatacrpdR1{t},2);
    end
    dif.length{i}=imresize(dif.length{i},[MeanBacLifed 1],'bilinear');
end

Tus.length=dif.length;



%% distance between tus spots and dif loci
%make a COM for dif
for j=1:size(dif.x,2)
    Xddd(:,j)=dif.x{j};
    Iddd(:,j)=dif.I{j};
end
for i=1:size(Xddd,1)
    dif.COM(i,1)=(Xddd(i,:)*Iddd(i,:)')/sum(Iddd(i,:));
end

% provide a distance from each Tus spot to the dif COM.
for j=1:size(Tus.x,2)
Tus.distance{j}=abs(nonzeros(Tus.x{j})-dif.COM(j,1));
end

%or do it with mean distance Tus (weighted)
for j=1:size(Tus.x,2)
    XTTT(:,j)=Tus.x{j};
    ITTT(:,j)=Tus.I{j};
end
for i=1:size(XTTT,1)
    Tus.COM(i,1)=(XTTT(i,:)*ITTT(i,:)')/sum(ITTT(i,:));
    Tus.xmean(i,1)=mean(nonzeros(XTTT(i,:)));
    Tus.xmeanstd(i,1)=std(nonzeros(XTTT(i,:)));
end



toc
