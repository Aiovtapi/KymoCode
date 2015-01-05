function dum=Processing00_TwoDSpot_ImageAnalyzer(exp)
%Two-dim analysis of EcoliData for bacterial analysis; Jacob Kerssemakers,
%TNW-BN-ND lab 2012; developed for Charl Moolman

if nargin<1, exp='TEST';end

actions.loaddatabase=1; %default=1 (new analysis)

initval=A001_Images_Set_Experiment(exp);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------


%---------------------------------------------

%work through all replication cycles------------------------
[~,chan_no]=size(S);
difs=[]; count=0;
for i=1:chan_no  %for each channel
Rep=S(i).channels.ReplicationCluster;
ManRep=M(i).channels.RepClicks;
chanstk_FL=S(i).channels.chanstk_FL;
kymo_FL=S(i).channels.kymo_FL;
[~,repno]=size(Rep);

for j=1:repno; %for each bacterium
'chan_no to go',  chan_no-i
'repno to go' ,repno-j
ThisRep=Rep(j);         %Current replication
ThisManRep=ManRep(j);

if strcmp(ThisManRep.fate, 'disassembled')
 
%fetch extra 'extended' data just before and right after the manually determined
%initiation and termination time points
ThisRep=Processing_Get_InitTerm_Extensions(ThisRep,initval);

%Do extensive 1D and 2D spot analysis. result: %pre-fit, final fit: 
%[X0,X1,Y0,Y1, Background amplitude,Peak0, Peak1, spotno]
[prefits,finalfits,spotnos,nofits,spotexcess]=Processing_ClusterLife(ThisRep,initval,chanstk_FL); 
%--------------------------------------------------------------------------

%Update the database with analysis result
%Add extended position data
S(i).channels.ReplicationCluster(j).PosKyTracCom.frames_ext=ThisRep.PosKyTracCom.frames_ext;  
S(i).channels.ReplicationCluster(j).PosKyTracCom.trackpos_ext=ThisRep.PosKyTracCom.frames_ext;  

%Fluorecencent spot: Add area summation data
S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontent_roi_ext=spotexcess;
S(i).channels.ReplicationCluster(j).FluoPropsGen.sumspot1_ext=nofits(:,1);
S(i).channels.ReplicationCluster(j).FluoPropsGen.sumspot2_ext=nofits(:,2);
S(i).channels.ReplicationCluster(j).FluoPropsGen.sumspotall_ext=nofits(:,3);

%Fluorecencent spots: add 2x 1D Gauss fit (along X) plus Y estimates
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.X0=prefits(:,1); %X0
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.X1=prefits(:,2); %X1
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.Y0=prefits(:,3); %Y0
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.Y1=prefits(:,4); %Y1
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.Bck=prefits(:,5);%Background
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.Pk1=prefits(:,6); %pixel contents spot 1
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.Pk2=prefits(:,7); %pixel contents  spot 2
S(i).channels.ReplicationCluster(j).Pos2DPreTrac.spots=spotnos; %#number of spots 

%Fluorecencent spots: add 2D Gauss fit 
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.X0=finalfits(:,1); %X0
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.X1=finalfits(:,2); %X1
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.Y0=finalfits(:,3); %Y0
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.Y1=finalfits(:,4); %Y1
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.Bck=finalfits(:,5);%Background
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.Pk1=finalfits(:,6); %pixel contents  spot 1
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.Pk2=finalfits(:,7); %pixel contents  spot 2
S(i).channels.ReplicationCluster(j).Pos2DFinTrac.spots=spotnos; %#number of spots      
%--------------------------------------------------------------------------

%plot menu for diagnosis---------------------------------------------------
X=ThisRep.PosKyTracCom.frames_ext;
%Y0=nofits(:,3);                             %spot count, areafit
Y2=(prefits(:,6)+prefits(:,5));            %spot count 1D fit
Y3=(finalfits(:,6)+finalfits(:,5))/pi;      %spot count 2D fits
%Y0=S(i).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontent_ext;
close(gcf);
scrsz = get(0,'ScreenSize');
figure('Position',[100 100 scrsz(3)/1.3 scrsz(4)/1.3]);

if S(i).channels.RepClicks(j).name>0  %no starters

mark=[-0.2*max(Y2) 1.2*max(Y2)];
fr0=S(i).channels.RepClicks(j).PosClick.firstframe;
fr1=S(i).channels.RepClicks(j).PosClick.lastframe;

subplot(1,2,1);
plot(X,Y2,'-o'); hold on
plot(X,Y3,'ko-');

plot([fr0 fr0],mark, 'k-');
plot([fr1 fr1],mark, 'k-');
axis([X(1) X(end) mark]);
title('total intensity  and spot(s) intensity, 1D fits')
xlabel('time, frames');
ylabel('intensity, pixel counts');
%legend('total intensity' ,'spots via summed backbone','clicked times');

% subplot(1,2,2);
% dum=Processing_Map_Fluorescence(ThisBac,ThisRep,chanstk_FL,kymo_FL,initval);

pause(0.5); %[~]=ginput(1);
%------------------------------------------------------
end
end
end
end
% 
save(outname, 'S', '-append');
disp('done');


