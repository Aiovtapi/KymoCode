function A060_RepliCluster00_Click(user,exp)

if nargin<2
    exp='Exp001_DnaN_TUS_dif_01092016_M';
end
if nargin<1
    user = 'MarkPC';
end

%Analyze_ReplicationCluster
%ReplicationCluster
%-------------------------------------------------------------------------
tic
close all;

% exp='001_DnaN_TUS_dif_30122014_difsignal';

actions.init=1;        %default=1 (new analysis)
actions.clickcycles=1; %default=1 (new analysis) else load
actions.divisionoverlay=0;
actions.eliminate_growth=1;

peak_tol=0;   %sensitivity for detecting peaks

initval=A001_Images_Set_Experiment(user,exp);

chans=initval.channelno;
for ch=1:chans
close all
display(strcat('chans to go = ',num2str(chans-ch)));
DnaNIdx=find(ismember(initval.viewchan,initval.DnaNchan));
Channelpath=char(strcat(initval.basepath,initval.nms{ch}{DnaNIdx},'.mat'));
load(Channelpath, 'chanstk_BF','chanstk_FL','endpoints', 'kymo_BF','kymo_FL','presets');
[r,c]=size(kymo_BF); 

if actions.divisionoverlay
BW=MakeBinaryEdgeImagefrom(kymo_BF,peak_tol);
kymo_FL_click(:,:,1)=kymo_FL*255/range(kymo_FL(:))+0.5*BW*255;
else
kymo_FL_click=kymo_FL;
end

%static_kymo=Processing_Straighten_Growth(kymo_FL,initval);  %TEST

kymoprops.width=c;
kymoprops.duration=r;
kymoprops.zoom=initval.zoom;  %used for clicking

if actions.init
    titl=strcat('Click Start positions, then right-click channel end position');
    
    [ReplicationCluster,RepClicks]=RepliCluster_Init(ch,kymo_FL_click,kymoprops, initval,actions);  %Init start spots
end
if actions.clickcycles
   figure(1)
[ReplicationCluster,RepClicks]=RepliCluster_GetCyclePoints(ch,ReplicationCluster,RepClicks,kymoprops, kymo_FL_click,initval,actions);
    %Detect (manually) start and end points of a replication spot 'cluster';
    %that may consist of one or two spots but is considered to represent one
    %replication cycle
    
else
%load  database--------------------------------------------------
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
load(outname_usr,'M');
[RepClicks,ReplicationCluster]=LoadCyclePoints(M,ch,initval);
end

ChanNum=size(initval.viewchan,2);

for i=1:ChanNum
kymoprops.WorkspaceOutName=char(initval.nms{ch}{i}); 
outname=strcat(initval.basepath,kymoprops.WorkspaceOutName);
save(outname, 'initval', 'RepClicks', 'ReplicationCluster',  '-append');
end
end
toc


function [RepClicks,ReplicationCluster]=LoadCyclePoints(M,ch,initval);
%Get cycle points and some first estimates from stored data inestead of
%clicking;
RepClicks=M(ch).channels.RepClicks;
[~,reps]=size(RepClicks);
for rp=1:reps 
rp
ystart=RepClicks(rp).PosClick.firstframe;
yend=RepClicks(rp).PosClick.lastframe;
xstart=RepClicks(rp).PosClick.firstpos;
xend=RepClicks(rp).PosClick.lastpos;
poss=linspace(xstart,xend,yend-ystart+1);
frs=[ystart:1:yend];
ReplicationCluster(rp).PosKyTracCom.frames=frs;
ReplicationCluster(rp).PosKyTracCom.clickpos=poss;
end


function BW=MakeBinaryEdgeImagefrom(KBF,tol);
% This function transforms an image into an edge-detected, binary image;
% this is useful if peaks with strongly varying intensities exist. JK13
    
BW=0*KBF;
[r,c]=size(KBF);
for i=1:r
    prf=KBF(i,:);
    mxi=F110_Get1DpeaksFlatBottom(prf,tol); %get peaks per line
    BW(i,mxi)=1;
    %BW(i,2)=1;   %just to have points in every 'frame'
    %BW(i,end-2)=1;
    prf=BW(i,:);
    lm=length(mxi); 
end
    %Series of binary operations
  BW=bwmorph(BW,'clean');      
  BW=bwmorph(BW,'dilate',1);
  %BW=bwmorph(BW,'clean');
  %BW=bwmorph(BW,'erode',1);
  BW=bwmorph(BW,'skel', Inf); 

     P_Color(1-BW,500,1000, 'grey');
% [~]=ginput(1);



