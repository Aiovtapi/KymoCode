%B0060_GetDistributionOfSpotRatio
%Load database, 
%Collect spots fit result
%Per spot pair, store the ratio dimspot/brightspot (per definiton 0..1)
%bin the results
%JacobKers 2014----------------------------------

close all
clear all
%select the experiment database-----------------------------

%exp='TEST';
exp='001_DnaN_TUS_dif_30122014_DnaNsignal_test';
%exp='CM_DnaN-Dif-Gamma-ve-Position1_Series1';
titlelabel=exp;
actions.savedata=1;             %default=1 (new analysis)
actions.loaddatabase=1;         %default=1 (new analysis)

binz=35;


initval=A001_Images_Set_Experiment(exp);
outname=strcat(initval.basepath,initval.outname);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------

[~,chan_no]=size(S);
[good, bad]=Processing_Measure_GoodDatabase(S)
totalbaccount=0;

ratio_ax=linspace(0,1,binz);
%arrays for collecting averaged data
%-------------------------------------

for ch=1:chan_no  %for each channel
chan_no-ch;
Div=S(ch).channels.AutoDivision;
Rep=S(ch).channels.ReplicationCluster;

kymo_FL=S(ch).channels.kymo_FL;
stripmov_FL=S(ch).channels.chanstk_FL;
[~,bacno]=size(Div);
allratios=[];
alldistances=[];
for j=1:bacno  %for each bacterium 
% bacno-bc
bacname=S(ch).channels.RepClicks(j).name;

ThisRep=Rep(j);         %Current replication
ThisBac=Div(j);
frs= ThisBac.edges.frs; 
ThisBac.family.me
j

%check number of conditions
 ok1=strcmp(ThisBac.birthtype, 'OK');    %birth ok
 ok2= strcmp(ThisBac.divtype, 'OK') ;   ;%division ok
 ok3=ThisBac.edges.edgesok;  %edges ok
 ok4=(length(frs)>0) ;
 ok5=strcmp(S(ch).channels.RepClicks(j).fate, 'disassembled');
 ok6=S(ch).channels.AutoDivision(j).accepted;  %manual accept/reject
 ok7=isfield(ThisRep,'Pos2DPreTrac');
if ok1&ok2&ok3&ok4&ok5&ok6&ok7
    totalbaccount=totalbaccount+1;
    Lspot=length(ThisRep.Pos2DPreTrac.contentspot1);
    bigspot=zeros(Lspot,1);
    littlespot=zeros(Lspot,1);
    thisbacdistances=zeros(Lspot,1);
    %Collect ratios for this bacterium
    for sp=1:Lspot
        switch 2
            case 1
            bigspot(sp)=max([ThisRep.Pos2DPreTrac.contentspot1(sp),ThisRep.Pos2DPreTrac.contentspot2(sp)]);
            %
            littlespot(sp)=min([ThisRep.Pos2DPreTrac.contentspot1(sp),ThisRep.Pos2DPreTrac.contentspot2(sp)]);
            thisbacdistances(sp) = abs(ThisRep.Pos2DPreTrac.X0(sp)-ThisRep.Pos2DPreTrac.X1(sp));
            case 2
            bigspot(sp)=max([ThisRep.Pos2DFinTrac.contentspot1(sp),ThisRep.Pos2DFinTrac.contentspot2(sp)]);
            littlespot(sp)=min([ThisRep.Pos2DFinTrac.contentspot1(sp),ThisRep.Pos2DFinTrac.contentspot2(sp)]);
            thisbacdistances(sp) = abs(ThisRep.Pos2DFinTrac.X0(sp)-ThisRep.Pos2DFinTrac.X1(sp));
        end
    end

    thisbacratios=littlespot./bigspot;      %ratios for this bacterium
    allratios=[allratios ; thisbacratios];  %collect all
    alldistances=[alldistances ; thisbacdistances];  %collect all
end
end
end


%-------------------------------------------------------------------------
%Here, we make a histogram of the spot ratio
%result.
%-------------------------
spothist=hist(allratios,ratio_ax);

figure; 
bar(ratio_ax,spothist,1,'r'); hold on;
title(strcat('Spot Ratios, ','N=',num2str(totalbaccount)));
xlabel('smallspot/bigspot, [-]');
ylabel('counts');
%-------------------------------------------------

figure;
length(allratios)
length(alldistances)
plot(allratios,alldistances,'o');

