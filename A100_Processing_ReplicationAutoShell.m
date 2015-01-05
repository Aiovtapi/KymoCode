%Processing_ReplicationAutoShell
%This script collects a series of automatic analysis steps
%JacobKers 2013----------------------------------------------------

tic
exp='002_DnaN_TUS_dif_21112014_DnaNsignal';
secondexp='002_DnaN_TUS_dif_21112014_TUSsignal';

%%-------------------------------------
%First, perform Center-off mass tracking on clusters starting at time
%points indicated by users 
disp('cleaning, quick tracking...');
if 1, RepliCluster00_TrackandcleanQuick(exp); end

%%-----------------------------------------------------------
%Next, Collect all channel data in one big database (just an administrative
%step)
disp('collecting...');
if 1, Processing_Collect_DataBases(exp); end

%%------------------------------------------------------------
%In this step, the moments of birth and division are deteced (from the brightfield data)  associated 
%with a replication cycle. In between these, the edges are detected time-point wise; 
%these points are cleaned from erroneous detections and used for
%(time-position) fits on the positions of this bacterium's edges 
disp('find division times...');
if 1, Processing_Find_Division_Times(exp); end %NB:still need to put off ginput for BW traces

%--------------------------------------------------------------------------
%Next, Get various fluorescence props like total fluorescence count, median excess
%count (a robust spots count estimate )
disp('adding general fluorescence info');
if 1, Processing_AnalyzeDivReptimingAuto(exp); end  

%--------------------------------------------------------------------------
%Next, a more detailed analysis on tthe precise times of initiation and
%termination (as opposed to the manual clicks) based on step fittng the
%spot focii signal
disp('finding init and ter......');
if 1, Processing_InitTer_analysisAuto(exp); end

%--------------------------------------------------------------------------
%Now, a detailed (and time consuming) analysis on the individual focci, based on first 1D-double
%Gaussian fitting, then a full double 2D Gaussian fit.

if 1, 
    disp('adding spot fluorescence info');
    Processing00_TwoDSpot_ImageAnalyzerAuto(exp); 
end

if 1, 
    disp('adding spot fluorescence info');
   Processing00_TwoDSpot_ImageAnalyzerAuto_SecondLabel(exp,secondexp); 
end
toc
%--------------------------------------------------------------------
%Next, a manual evaluation of bacterium cycles, based on fluorescence
%position-time graphs. This part includes automatic cleaning steps (which
%is preferred) before the user applies final judgment.
if 1, A111_Processing_ManualAcceptRejectBacteriaAutoDiv(exp); end



