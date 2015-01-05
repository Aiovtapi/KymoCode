function Processing_AnalyzeDivReptiming(exp)
%Get some generalfluorescent properties of the fluorescence signal
%JacobKers 2012

actions.plot=0;


close all
actions.savedata=1;         %default=1 (new analysis)
actions.loaddatabase=1;     %default=1 (new analysis)
actions.plot=0;             %default=0 (quick run)

if nargin<1, 
    %exp='TEST';
    %exp='CM_HFR1_JK' ; 
    exp='2a';
end

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

goodcount=0;
[~,chan_no]=size(S);
%chan_no=2;
difs=[]; count=0;
for i=1:chan_no  %for each channel
chan_no-i
Div=S(i).channels.Division;
Rep=S(i).channels.ReplicationCluster;
kymo_FL=S(i).channels.kymo_FL;
stripmov_FL=S(i).channels.chanstk_FL;
[~,bacno]=size(Div);
[~,repno]=size(Rep);
for j=1:bacno  %for each bacterium
'bacno', bacno-j
spotno=[];
ThisBac=Div(j);

linkedrep=ThisBac.linkedrep;

if linkedrep>0, ThisRep=Rep(linkedrep);end

frs=ThisBac.PosKyTrac.frames; 
divtime=frs(end); 
birthtime=frs(1);     
lft=ThisBac.PosKyTrac.left;
rht=ThisBac.PosKyTrac.right;

if linkedrep>0
count=count+1;
ThisRep=Rep(linkedrep);
repfrs=ThisRep.PosKyTracCom.frames;
reppos=ThisRep.PosKyTracCom.trackpos;
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

%obtain averaged profile
prof=kymo_FL(fri,idx); prof=prof-min(prof)+0.001;
lixx=length(idx); shft=(hor-lixx)/2;
idx2=round([1:lixx]+shft); %centered indices;
pic(t,idx2)=prof; 


%get fluorescent  area ful bacterium----------------
flipstrip=fliplr(stripmov_FL(:,:,fri));
FL=flipstrip(:,idx);
[r,c]=size(FL);

%get fluorescent  area Spot ROI----------------
FL_roi=flipstrip(:,idx_roi);
[r_roi,c_roi]=size(FL_roi);

%do image analysis, no fitting on full bacterial area--------------------
fluo=Fluorescence_PatternAnalysis(FL); 
fluo_roi=Fluorescence_PatternAnalysis(FL_roi); 


%transfer data
%...from full bacterium area
fluodarklevel(t)=fluo.darklevel;
fluosignalcontent(t)=fluo.signalcontent;
fluospinelevel(t)=fluo.maxspinelevel; %this is later used as bacterium background
fluospotscontent(t)=fluo.spotscontent;
peakval(t)=fluo.peakval;


%---from limited region (to avoid offspring spots)


if c>1,peakx(t)=fluo.peakx+idx(1)-1;else peakx(t)=1;end
peaky(t)=fluo.peaky;
end

%define a 'compact' slot between the originally clicked points-------
lf=length(repfrs);
fluodarklevel_compact=[];
fluosignalcontent_compact=[];
fluospinelevel_compact=[];
fluospotscontent_compact=[];
fluopeakval_compact=[];
fluopeakx_compact=[];
fluopeaky_compact=[];
for ti=1:lf
idx=find(repfrs(ti)==frs_ext); 
fluodarklevel_compact(ti)=fluodarklevel(idx);
fluosignalcontent_compact(ti)=fluosignalcontent(idx);
fluospinelevel_compact(ti)=fluospinelevel(idx);
fluospotscontent_compact(ti)=fluospotscontent(idx);
fluopeakval_compact(ti)=peakval(idx);
fluopeakx_compact(ti)=peakx(idx);
fluopeaky_compact(ti)=peaky(idx);
end

%Store the  fluorescence properties in the
%database-------------
%S(i).channels.ReplicationCluster(linkedrep).ThisRep.PosKyTracCom.frames

S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.darklevel=fluodarklevel_compact;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.darklevel_ext=fluodarklevel;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalcontent=fluosignalcontent_compact;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalcontent_ext=fluosignalcontent;

S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalpeakval=fluopeakval_compact;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalpeakval_ext=peakval;

S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.fluospinelevel=fluospinelevel_compact;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.fluospinelevel_ext=fluospinelevel;

S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.fluospotscontent=fluospotscontent_compact;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.fluospotscontent_ext=fluospotscontent;

S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalpeakx=fluopeakx_compact;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalpeakx_ext=peakx;

S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalpeaky=fluopeaky_compact;
S(i).channels.ReplicationCluster(linkedrep).FluoPropsGen.signalpeaky_ext=peaky;
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
plot(rht-middiv+hor/2+1.5,frs-frs_ext(1)+1,'-o'); hold on
plot(lft-middiv+hor/2+1.5,frs-frs_ext(1)+1,'-o');

title(strcat('Bacterium',num2str(ThisBac.name)));
ylabel('time, frames');
xlabel('position, pixels');
axis([1 hor 1 ver]);

fluo2=fluosignalcontent;

subplot(1,2,2);
plot(fluo2(1:ver),frs_ext, '-o');
axis([1 1.3*max(fluo2) frs_ext(1) frs_ext(end)]);
title('Fluorescence')
xlabel('Total Fluorescence, a.u.');
ylabel('time, frames');


close(gcf);
end
%------------------
end
end
end

%Finally, save stuff
if actions.savedata
    save(outname,'S', '-append');
end
display('done');
end



function fluo=Fluorescence_PatternAnalysis(FL);
%This function returns some fluorescent properties, based on simple
%assumption of the fluorescent pattern of the bacterium
[r,c]=size(FL);
if c>1
%estimate summed background intensity over image

fluo.darklevel=(mean(mean(FL(1:2,:)))+mean(mean(FL(r-1:r,:))))/2;
%'local background' is defined as the average of the outer two image lines'
fluosumcurve=sum(FL)-r*fluo.darklevel;
fluomaxcurve=max(FL)-fluo.darklevel;

fluo.signalcontent=sum(fluosumcurve);
fluo.sumspinelevel=median(fluosumcurve);  %total represents total counts
fluo.maxspinelevel=median(fluomaxcurve);  %level represents backbone level bacterium


%Get rough estimate of most prominent features
[fluomaxval,xidx]=max(fluomaxcurve);
[~,yidx]=max(max(FL'));


sel=find(fluosumcurve>fluo.sumspinelevel); 
fluo.spotscontent=sum(fluosumcurve(sel)-fluo.sumspinelevel);

fluo.peakval=fluomaxval-fluo.maxspinelevel;  %estimate peak-only value
fluo.peaky=yidx;
fluo.peakx=xidx;  %xcoordinate in absolute indices kymoline
else
fluo.darklevel=0;
fluo.signalcontent=0;
fluo.maxspinelevel=0;
fluo.sumspinelevel=0;
fluo.spotscontent=0;
fluo.peakval=0;
fluo.peaky=1;
fluo.peakx=1;
end

if 0
fluo
plot(fluosumcurve, 'k-o');hold on; 
dummy=0*fluosumcurve;
%plot(dummy+fluo.darklevel,'-');
plot(dummy+fluo.sumspinelevel,'r-');

axis([0 length(sum(FL)) 0 1.2*max(fluosumcurve)]);

spotpercent=fluo.spotscontent/fluo.signalcontent*100
[~]=ginput(1);
close(gcf);
end

dum=1;


