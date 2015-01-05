function dum=A111_Processing_ManualAcceptRejectBacteriaAutoDiv(exp)
%This is a simple loop fetching properties of associated division and
%replication cycles, allowing the user to accept or reject the result

%Since this is hand work, it is also saved separately in an user-input 'M' database, for use on different sets.
%JacobKers 2012

close all;
actions.dosave=1;  %default=1 (new analysis) Use this when props need to be stored   
testNur=1
%select the experiment database-----------------------------
if nargin<1
exp='B_DnaN_TUS_dif_21112014_DnaNSignal';
%exp='2a';
%exp='CM_HFR1_JK'
%exp='DnaN_20msExp_10minFramRate Series 1'
%exp='20121214_DnaN_20msExpTIme_TIffs_004';
%exp='20130221_HigherFrameRateMeasurement';
%exp='dDnaX-mYpet-MM' ;
end

initval=A001_Images_Set_Experiment(exp);
outnameS=strcat(initval.basepath,initval.outname);
outnameM=strcat(initval.basepath,initval.outname_usr);
load(outnameS,'S');
load(outnameM,'M');
%--------------------------------------------------
ThisIsSizeS = size(S)

[~,chan_no]=size(S)
%[chan_no,~]=size(S)


%First, quick count of database
Nbac=Processing_Measure_Database(S)

%Loop through all bacteria-----------------------------
count=0;
count.good=0;
count.bad=0;
for i=1:chan_no  %for each channel
chan_no-i

Div=S(i).channels.AutoDivision;
Rep=S(i).channels.ReplicationCluster;
RepC=S(i).channels.RepClicks;

kymo_FL=S(i).channels.kymo_FL;
stripmov_FL=S(i).channels.chanstk_FL;
[~,bacno]=size(Div);
for j=1:bacno  %for each bacterium   
display('Bacteria to go'), Nbac-count.good-count.bad
ThisBac=Div(j); 
ThisRep=Rep(j);

    ok1=strcmp(ThisBac.birthtype, 'OK');    %birth ok
    ok2= strcmp(ThisBac.divtype, 'OK');    ;%division ok
    ok3=ThisBac.edges.edgesok;  %edges ok   

    if ok2&ok3       
        dum=Processing_Map_FluorescenceAutoDiv(ThisBac,ThisRep,stripmov_FL,kymo_FL,initval);
        title('Left-click if OK; Right-click if reject');
        [~,~,but]=ginput(1);   
    else
        but=3;
    end
    
FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'ManualAcceptReject/Accept/'));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'ManualAcceptReject/Accept/'));
end

FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'ManualAcceptReject/Reject/'));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'ManualAcceptReject/Reject/'));
end
    

switch but
case 1
S(i).channels.AutoDivision(j).accepted=1;
M(i).channels.RepClicks(j).accepted=1;
count.good=count.good+1;
FigureExistence = findobj('type','figure');
if FigureExistence ==1
    h=gcf;
    print(h, '-dpng', '-r300',strcat(initval.basepath,initval.FiguresFolder,'ManualAcceptReject/Accept/','Ch',int2str(i),'_BacNo',int2str(j)));
end;
%'accept'
case 3
S(i).channels.AutoDivision(j).accepted=0;
M(i).channels.RepClicks(j).accepted=0;
count.bad=count.bad+1;
FigureExistence = findobj('type','figure');
if FigureExistence ==1
    h=gcf;
    print(h, '-dpng', '-r300',strcat(initval.basepath,initval.FiguresFolder,'ManualAcceptReject/Reject/','Ch',int2str(i),'_BacNo',int2str(j)));
end;
end
close(gcf);
end
but
end

if actions.dosave

save(outnameS, 'S','M', '-append');
save(outnameM, 'M', '-append');
disp('done');
end
dum=1;

