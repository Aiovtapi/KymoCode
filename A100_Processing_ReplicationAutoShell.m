%Processing_ReplicationAutoShell
%This script collects a series of automatic analysis steps
%JacobKers 2013----------------------------------------------------

function A100_Processing_ReplicationAutoShell(user,exp,WorkspaceOutname)

if nargin<2
    exp='Exp001_DnaN_TUS_dif_01092016_M';
end
if nargin<1
    user = 'MarkPC';
end

tic

initval=A001_Images_Set_Experiment(user,exp);

if nargin == 3;
    initval.nms = WorkspaceOutname;
end

DnaNIdx=find(ismember(initval.viewchan,initval.DnaNchan));
ChanNum=size(initval.viewchan,2);
Dummy=linspace(1,ChanNum,ChanNum);
Dummy(Dummy==DnaNIdx)=[];


%%-------------------------------------
%First, perform Center-off mass tracking on clusters starting at time
%points indicated by users 
disp('cleaning, quick tracking...');
if 1, RepliCluster00_TrackandcleanQuick(exp,user,DnaNIdx,initval.nms); end

%%-----------------------------------------------------------
%Next, Collect all channel data in one big database (just an administrative
%step)


disp('collecting...');
if 1, Processing_Collect_DataBases(exp,user,DnaNIdx,DnaNIdx,initval.nms); end

%%------------------------------------------------------------
%In this step, the moments of birth and division are detected (from the brightfield data)  associated 
%with a replication cycle. In between these, the edges are detected time-point wise; 
%these points are cleaned from erroneous detections and used for
%(time-position) fits on the positions of this bacterium's edges 
disp('find division times...');
if 1, Processing_Find_Division_Times(exp,user,DnaNIdx); end %NB:still need to put off ginput for BW traces

%--------------------------------------------------------------------------
%Next, Get various fluorescence props like total fluorescence count, median excess
%count (a robust spots count estimate )
disp('adding general fluorescence info');
if 1, Processing_AnalyzeDivReptimingAuto(exp,user,DnaNIdx); end  

%--------------------------------------------------------------------------
%Next, a more detailed analysis on tthe precise times of initiation and
%termination (as opposed to the manual clicks) based on step fittng the
%spot focii signal
disp('finding init and ter......');
if 1, Processing_InitTer_analysisAuto(exp,user,DnaNIdx); end

%--------------------------------------------------------------------------
%Now, a detailed (and time consuming) analysis on the individual foci, based on first 1D-double
%Gaussian fitting, then a full double 2D Gaussian fit.

if 1, 
    disp('adding spot fluorescence info');
    Processing00_TwoDSpot_ImageAnalyzerAuto(exp,user,DnaNIdx,DnaNIdx); 
end



% This can be faster as only Processing_Clusterlife is needed for the other
% channels. Gaussian fitting is now done 3 times. 
for N=Dummy;
    Processing_Collect_DataBases(exp,user,N,DnaNIdx,initval.nms);
    Processing00_TwoDSpot_ImageAnalyzerAuto(exp,user,N,DnaNIdx);
end

% The folders of the bacpics will be renamed so that the indices of the
% folders corresponds to the foldernames. 
disp('Renaming Bacfolders')
if 1, Processing_Relabel_Bacfolders(exp,user); end

toc

%--------------------------------------------------------------------
%Next, a manual evaluation of bacterium cycles, based on fluorescence
%position-time graphs. This part includes automatic cleaning steps (which
%is preferred) before the user applies final judgement.
%if 1, A111_Processing_ManualAcceptRejectBacteriaAutoDiv(exp); end



