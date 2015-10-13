%Processing_ReplicationAutoShell
%This script collects a series of automatic analysis steps
%JacobKers 2013----------------------------------------------------

tic
exp='001_DnaN_TUS_dif_30122014_DnaNsignal';
expFLchan2='001_DnaN_TUS_dif_30122014_difsignal';

initval=A001_Images_Set_Experiment(exp);
%%-------------------------------------
%First, perform Center-off mass tracking on clusters starting at time
%points indicated by users 

if 1, 
    disp('cleaning, quick tracking...');
    RepliCluster00_TrackandcleanQuick(exp); 
end

%%-----------------------------------------------------------
%Next, Collect all channel data in one big database (just an administrative
%step)

if 1, 
    disp('collecting...');
    Processing_Collect_DataBases_2FLchan(exp,expFLchan2); 
end

%%------------------------------------------------------------
%In this step, the moments of birth and division are deteced (from the brightfield data)  associated 
%with a replication cycle. In between these, the edges are detected time-point wise; 
%these points are cleaned from erroneous detections and used for
%(time-position) fits on the positions of this bacterium's edges 

if 1, 
    disp('find division times...');
    Processing_Find_Division_Times(exp); 
end %NB:still need to put off ginput for BW traces

%--------------------------------------------------------------------------
%Next, Get various fluorescence props like total fluorescence count, median excess
%count (a robust spots count estimate )

if 1, 
    disp('adding general fluorescence info');
    Processing_AnalyzeDivReptimingAuto(exp); 
end  

%--------------------------------------------------------------------------
%Next, a more detailed analysis on tthe precise times of initiation and
%termination (as opposed to the manual clicks) based on step fittng the
%spot focii signal

if 1, 
    disp('finding init and ter......');
    Processing_InitTer_analysisAuto(exp); 
end

%--------------------------------------------------------------------------
%Now, a detailed (and time consuming) analysis on the individual focci, based on first 1D-double
%Gaussian fitting, then a full double 2D Gaussian fit.
%--------------------------------------------------------------------
%Next, a manual evaluation of bacterium cycles, based on fluorescence
%position-time graphs. This part includes automatic cleaning steps (which
%is preferred) before the user applies final judgment
if 1, 
    disp('waiting for your verdicts....');
    A111_Processing_ManualAcceptRejectBacteriaAutoDiv(exp); 
end

if 1, 
    disp('adding spot fluorescence info; pass 1.....');
    Processing00_TwoDSpot_ImageAnalyzerAuto_FL2chan(exp,1); %DnaN channel
    disp('adding spot fluorescence info; pass 2.....');
    Processing00_TwoDSpot_ImageAnalyzerAuto_FL2chan(exp,2); %dif channel
end






