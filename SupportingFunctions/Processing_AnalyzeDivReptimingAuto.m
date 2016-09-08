function Processing_AnalyzeDivReptimingAuto(exp,user,ColourIdx)
%Get some generalfluorescent properties of the fluorescence signal
%JacobKers 2012


close all


if nargin<1, 
    exp='TEST';
    %exp='CM_HFR1_JK' ; 
    %exp='CM_DnaN_37Deg_Series1002';
end

initval=A001_Images_Set_Experiment(user,exp);


actions.savedata=1;         %default=1 (new analysis)
actions.loaddatabase=1;     %default=1 (new analysis)
actions.plot=initval.plotintermediateresults;             %default=0 (quick run)


%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname{ColourIdx});
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------

goodcount=0;
%[~,chan_no]=size(S);
ThisIsSizeS = size(S);
%[chan_no,~]=size(S)%%I THINK THIS WAS AN ERROR
[~,chan_no]=size(S)
%chan_no=2;
difs=[]; count=0;
for i=1:chan_no  %for each channel
chan_no-i
Div=S(i).channels.AutoDivision;
Rep=S(i).channels.ReplicationCluster;
kymo_FL=S(i).channels.kymo_FL;
stripmov_FL=S(i).channels.chanstk_FL;
[~,bacno]=size(Div);
for j=1:bacno  %for each bacterium
'bacno', bacno-j
spotno=[];
ThisBac=Div(j);
ThisRep=Rep(j);
frs= ThisBac.edges.frs; 

%check number of conditions
 ok1=strcmp(ThisBac.birthtype, 'OK');    %birth ok
 ok2= strcmp(ThisBac.divtype, 'OK');    ;%division ok
 ok3=ThisBac.edges.edgesok;  %edges ok
 ok4=(length(frs)>0);
 %ok5=ThisBac.accepted;
 
%  indivbac = j
% chann=i
% %repno=repno
% bacno=bacno
%  brithOk = ok2
%  divok = ok3
%  lengthok = ok4
% terok = ok6
 
 
if ok2&ok3&ok4
divtime=frs(end); 
birthtime=frs(1);     

rht= ThisBac.edges.leftfit+1;  % 
lft= ThisBac.edges.rightfit-1;
% rht= ThisBac.edges.leftfit;
% lft= ThisBac.edges.rightfit;


count=count+1;
repfrs= ThisBac.edges.repfrs;
reppos=ThisBac.edges.reppos;
%diffa=ThisRep.Pos2DPreTrac.X1-ThisRep.Pos2DPreTrac.X0;

%some more props
middiv=(rht+lft)/2;  %midline bactrium

%select a time slot encompassing replication and division cycle
%-----------------------------------------------------
xt=initval.extension;  %extension before and after replication and division times
hifr=repfrs(end)+xt;
lofr=repfrs(1)-xt;


frs_ext=[lofr:hifr];

%Fits to bacterium position------------------------------------------
%fits to midline bacterium, second order
ppM=polyfit(frs,middiv,2);
fit_mid=ppM(1)*(frs_ext).^2+ppM(2)*frs_ext+ppM(3); 

%fits RELATIVE to MIDLINE
ppL=polyfit(frs,lft-middiv,1);   
fit_lft=round(ppL(1)*(frs_ext)+ppL(2));  

ppR=polyfit(frs,rht-middiv,1);       
fit_rht=round(ppR(1)*(frs_ext)+ppR(2)); 
%-----------------------------------------------------------------

%using fit results, we now fetch a fluorescence profile from the kymograph:
%----------------------------------------------
hor=max(fit_rht-fit_lft)+4;
ver=length(frs_ext);
pic=zeros(ver,hor);
size(pic);
[kr,kc]=size(kymo_FL);
fluosignalcontent=[];
fluospinelevel=[];
fluospotscontent=[];
fluodarklevel=[];
peakval=[];
peakx=[];
peaky=[];
fluolevel_sumspine=[];

%for t=1:ver
for t=1:ver
fri=max(frs_ext(t),1); fri=min(fri,kr);  %inbound frame number

%Get boundary edges: bacterium edges
loix=max(fit_lft(t)+fit_mid(t), 1);  %inbound low pos index
hiix=min(fit_rht(t)+fit_mid(t), kc); %inbound hi pos index
idx=round([loix:hiix]);% indices for kymograph line

%Get boundary edges: spot region
loix_roi=max(fit_mid(t)-initval.spotRoiHW, 1);  %inbound low pos index
hiix_roi=min(fit_mid(t)+initval.spotRoiHW, kc); %inbound hi pos index
idx_roi=round([loix_roi:hiix_roi]);% indices for kymograph line




%get fluorescent  area ful bacterium----------------
flipstrip=fliplr(stripmov_FL(:,:,fri));
FL=flipstrip(:,idx);
[r,c]=size(FL);

%get fluorescent  area Spot ROI----------------
FL_roi=flipstrip(:,idx_roi);
[r_roi,c_roi]=size(FL_roi);


%do image analysis, no fitting on full bacterial area--------------------
fluo=Processing_Fluorescence_PatternAnalysis(FL); 
fluo_roi=Processing_Fluorescence_PatternAnalysis(FL_roi); 

%obtain maximized profile
%prof=kymo_FL(fri,idx); prof=prof-min(prof)+0.001;



%%- I COMMENTED THIS OUT DATE: 24 NOVEMBER 2013, SINCE IT GAVE AN ERROR -
%%SOMETIMES
prof=fluo.curve_medianofmax;
lixx=length(idx); shft=(hor-lixx)/2;
idx2=round([1:lixx]+shft); %centered indices;
pic(t,idx2)=prof; 





%transfer data
%...from full bacterium area
fluodarklevel(t)=fluo.level_dark;
fluosignalcontent(t)=fluo.content_signal;
fluospinelevel(t)=fluo.level_medianofmax; %this is later used as bacterium background
fluospotscontent(t)=fluo.content_spots;
fluolevel_sumspine(t)=fluo.level_medianofsum;
peakval(t)=fluo.level_peak;
end




%define a 'compact' slot between the originally clicked points-------
lf=length(repfrs);
fluodarklevel_compact=[];
fluosignalcontent_compact=[];
fluospinelevel_compact=[];
fluospotscontent_compact=[];
fluopeakval_compact=[];

for ti=1:lf
idx=find(repfrs(ti)==frs_ext);
idx;
ti;
fluodarklevel_compact(ti)=fluodarklevel(idx);
fluosignalcontent_compact(ti)=fluosignalcontent(idx);
fluospinelevel_compact(ti)=fluospinelevel(idx);
fluospotscontent_compact(ti)=fluospotscontent(idx);
fluopeakval_compact(ti)=peakval(idx);
end

%Store the  fluorescence properties in the
%database-------------
%S(i).channels.ReplicationCluster(j).ThisRep.PosKyTracCom.frames
S(i).channels.ReplicationCluster(j).FluoPropsGen.frs_ext=frs_ext;
S(i).channels.ReplicationCluster(j).FluoPropsGen.frs=repfrs;
S(i).channels.ReplicationCluster(j).FluoPropsGen.darklevel=fluodarklevel_compact;
S(i).channels.ReplicationCluster(j).FluoPropsGen.darklevel_ext=fluodarklevel;
S(i).channels.ReplicationCluster(j).FluoPropsGen.signalcontent=fluosignalcontent_compact;
S(i).channels.ReplicationCluster(j).FluoPropsGen.signalcontent_ext=fluosignalcontent;
S(i).channels.ReplicationCluster(j).FluoPropsGen.fluosumspinelevel_ext=fluolevel_sumspine;

S(i).channels.ReplicationCluster(j).FluoPropsGen.signalpeakval=fluopeakval_compact;
S(i).channels.ReplicationCluster(j).FluoPropsGen.signalpeakval_ext=peakval;

S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospinelevel=fluospinelevel_compact;
S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospinelevel_ext=fluospinelevel;

S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontent=fluospotscontent_compact;
S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontent_ext=fluospotscontent;
%-------------------------------------------


picmx=max(pic(:));
sumFL=(sum(pic'))';

%NOTE: kymographs store only maxima perpendicular to channel line (to
%enhance spot contrast for detection)
%total fluorescence is stored in channel strip movies



%plot menu; can also be used to judge data! ---------------------------
if actions.plot
close(gcf);
subplot(1,2,1);

pcolor(pic); colormap hot; shading flat;  hold on
plot(rht-middiv+hor/2+1.5,frs-frs_ext(1)+1,'-o','LineWidth',2,'MarkerSize',6); hold on
plot(lft-middiv+hor/2+1.5,frs-frs_ext(1)+1,'-o','LineWidth',2,'MarkerSize',6);
 set(gca, 'fontsize', 26, 'linewidth', 2, 'fontweight', 'bold','box','off');
axis xy;
set(gca,'YDir','reverse');
ylabel('Time (frames)', 'fontsize', 26, 'fontweight', 'bold')
xlabel('Position (px)', 'fontsize', 26, 'fontweight', 'bold')


axis([1 hor 1 ver]);

fluo2=fluospotscontent./fluosignalcontent*100;

subplot(1,2,2);
plot(fluo2(1:ver),frs_ext,'-bo','LineWidth',4,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','w',...
                'MarkerSize',6);
set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold','box','off');
set(gca,'TickLength',[0.02 0.02]);
axis([1 1.3*max(fluo2) frs_ext(1) frs_ext(end)]);
title('Fluorescence')
axis xy;
set(gca,'YDir','reverse');
xlabel('Tot Fluor in cell (a.u.)', 'fontsize', 26, 'fontweight', 'bold');
ylabel('Time (fr)', 'fontsize', 26, 'fontweight', 'bold');
pause(0.5)


FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'AnalyzeDivRepTimeAuto/'));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'AnalyzeDivRepTimeAuto/'));
end

h=gcf;
print(h, '-dpng', '-r300',strcat(initval.basepath,initval.FiguresFolder,'AnalyzeDivRepTimeAuto/','Ch',int2str(i),'_BacNo',int2str(j)));


close(gcf);
end
%------------------
else
    
S(i).channels.ReplicationCluster(j).FluoPropsGen.frs_ext=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.frs=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.darklevel=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.darklevel_ext=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.signalcontent=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.signalcontent_ext=[];


S(i).channels.ReplicationCluster(j).FluoPropsGen.signalpeakval=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.signalpeakval_ext=[];

S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospinelevel=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospinelevel_ext=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.fluosumspinelevel_ext=[];

S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontent=[];
S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontent_ext=[];
    
    
    
end






end
end
%Finally, save stuff
if actions.savedata
    save(outname,'S', '-append');
end
display('done');




