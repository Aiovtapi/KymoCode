function dum=B0020_Processing_PlotIntensityPerLengthafoExpTimeAuto
%This is a simple loop fetching properties of associated division and
%replication cycles, allowing the user to accept or reject the result

%JacobKers 2012

close all;
clear;

actions.savedata=1;
%fluotype='All per length';
fluotype='Peaks';
%fluotype='Cytoplasm per length';

extras.remark=fluotype;
binax=[0:1:300];
%select the experiment database-----------------------------

%exp='TEST';
%exp='2a';
%exp='CM_HFR1_JK'
%exp='DnaN_20msExp_10minFramRate Series 1'
%exp='20121214_DnaN_20msExpTIme_TIffs_004';
%%exp='20130221_HigherFrameRateMeasurement';
%exp='dDnaX-mYpet-MM' ;
%exp='20130529_DnaX_5minFrameRate';
%exp='CM_DnaN-Dif-Gamma-ve-Position1_Series1';
%exp='D005_TIFFs_DnaN_Dif_004';
%exp='CM_DnaN_20140103_D003_DnaN_DelTus';
exp='001_DnaN_TUS_dif_30122014_DnaNsignal_test';

initval=A001_Images_Set_Experiment(exp);
initval.expno='001_DnaN_TUS_dif_30122014_DnaNsignal_test';

outnameS=strcat(initval.basepath,initval.outname);
outnameM=strcat(initval.basepath,initval.outname_usr);
load(outnameS,'S');
load(outnameM,'M');
%--------------------------------------------------
[~,chan_no]=size(S); 

%First, quick count of database
Nbac=Processing_Measure_Database(S);
startnos=zeros(Nbac,1);
Xdata=[];
Ydata=[];

%Loop through all bacteria-----------------------------
count=0;
for i=1:chan_no  %for each channel
chan_no-i
Div=S(i).channels.AutoDivision;
Rep=S(i).channels.ReplicationCluster;
kymo_FL=S(i).channels.kymo_FL;
stripmov_FL=S(i).channels.chanstk_FL;
[~,bacno]=size(Div);
for j=1:bacno  %for each bacterium
'bacno', bacno-j
ThisBac=Div(j);
ThisRep=Rep(j);
%Get spot data---------------
frs= ThisBac.edges.frs; 


%check number of conditions-------------------------------------
 %ok1=strcmp(ThisBac.birthtype, 'OK');    %birth ok
 ok2= strcmp(ThisBac.divtype, 'OK');    ;%division ok
 ok3=ThisBac.edges.edgesok;  %edges ok
 ok4=(length(frs)>0); 
 ok5=strcmp(S(i).channels.RepClicks(j).fate, 'disassembled');


if ok2&ok3&ok4&ok5  %&acc  
count=count+1;
display('Bacteria to go'), Nbac-count

%--------------------------------------------------------------------------
%Select the division data; find corresponding fluorescence data and sum
%this; divide it by bacterium length and plot
[Timax,IntperLength]=Processing_Map_Div_Fluorescence(ThisBac,ThisRep,stripmov_FL,kymo_FL,initval,fluotype);
Xdata=[Xdata Timax];
Ydata=[Ydata IntperLength];
%--------------------------------------------------------------------------

end
end
dum=1;
end
dum=Binfiller(Xdata,Ydata,binax,initval,extras);

% %Saving menu---------------------------------------------------------------
% if actions.savedata
%     %1) histograms
%      %2) binned averages
%      avdata=[binaxnw' loall avall hiall lopeaks avpeaks hipeaks loperc avperc hiperc];
%     filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',exp,variant,remark,'_SpotsAndTotals.txt');
%      dlmwrite(filname,avdata, 'delimiter', '\t');
%      %2) individual points per bin
%      x=(binax2+scatax-minsperframe/2)';
%      x1=binax2'-minsperframe/2;
%      binpoints=[x1(:) x(:) bincollector_all_labels(:) bincollector_peak_labels(:) bincollector_perc_labels(:)];
%     filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',exp,variant,remark,'_BinPoints.txt');
%      dlmwrite(filname,binpoints, 'delimiter', '\t');
% end

function [Timax,IntperLength]=Processing_Map_Div_Fluorescence(ThisBac,ThisRep,stripmov_FL,kymo_FL,initval,fluotype);
%Function makes a map of fluorescence and division props
    
%frames and positions------------------------------------------------------
frs=ThisBac.edges.frs;      
lft=ThisBac.edges.rightfit;  %CHECK THIS BUG
rht=ThisBac.edges.leftfit;
repfrs=ThisBac.edges.repfrs;
%------------------------------------------------------

lfr=length(frs);
c=0;
for i=1:lfr;                        %for each time point of division cycle
    repfr_idx=find(repfrs==frs(i)); %find the corresponding rep signal
    if ~isempty(repfr_idx)          %if it exists
        c=c+1;
        BacLength=rht(i)-lft(i);
        switch fluotype
            case 'All per length', 
                SumInt=ThisRep.FluoPropsGen.signalcontent(repfr_idx);  %All!
                IntperLength(c)=SumInt/BacLength;   
            case 'Peaks', 
                SumInt=ThisRep.FluoPropsGen.fluospotscontent(repfr_idx);  %spots
                IntperLength(c)=SumInt;
            case 'Cytoplasm per length', 
                SumInt=(ThisRep.FluoPropsGen.signalcontent(repfr_idx)-ThisRep.FluoPropsGen.fluospotscontent(repfr_idx))/BacLength; 
                IntperLength(c)=SumInt;%cytoplasma per length
        end
        
        
        Timax(c)=frs(i);
        
        %  %for peak val!
    end
end
avT=mean(Timax);
avIperL=mean(IntperLength);


if 0
subplot(1,2,1);
plot(Timax,IntperLength, '*', 'MarkerSize', 3); hold on;
plot(avT,avIperL, 'ro', 'MarkerSize', 7); hold on;
axis([0 400 0 5E4]);
subplot(1,2,2);
plot(Timax,IntperLength, '-ko');
%[~]=ginput(1);
end
dum=1;

