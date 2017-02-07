clc
clear all

%%

TraceNumber=11;

init.OSslash = '/';

fprintf('Select Data Folder');
init.datapath = uigetdir(pwd,'Select Data Folder');
    
if init.datapath == 0;
    init.datapath = '/Users/rleeuw/Work/Data/170111_Tus-SMcal/gain100/11';
end

fprintf('\nSelect Utrack Results Folder');
init.utrackpath = uigetdir(init.datapath,'Select utrack Folder');

if init.utrackpath == 0;
    init.utrackpath = '/Users/rleeuw/Work/Data/170111_Tus-SMcal/gain100/11/utrackResults';
end

    init.datapath = strcat(init.datapath,init.OSslash);
    init.utrackpath = strcat(init.utrackpath,init.OSslash);
    
%% load the detection, tracking and motion analysis files

TP='TrackingPackage/';

utrack.Detection=load(strcat(init.utrackpath,TP,'GaussianMixtureModels/Channel_1_detection_result.mat'));
% utrack.Motion=load(strcat(init.utrackpath,TP,'MotionAnalysis/channel_1.mat'));
utrack.tracks=load(strcat(init.utrackpath,TP,'tracks/Channel_1_tracking_result.mat'));


%% Get the mean background intensity (and std) from all detections made per frame

Background=LionLoadBackground(utrack);

Nframes=size(Background,1);

for i=1:Nframes
Bkg.mean(i,1)=mean(Background{i});
Bkg.std(i,1)=std(Background{i});
end

clear Background

% create a cell containing matrices with  for
% each track, in this format: [X Y Z AMP Xstd Ystd Zstd Ampstd BkgMean Bkgstd]

[TrackCell,FrameArray,Bkg]=LionLoadTrackData(utrack,Bkg);


%% Collective Trace Plot

% all traces need to be as long as the experiment
Ntracks=size(TrackCell,1);
Trace=cell(Ntracks,1);
Frames=cell(Ntracks,1);


for i=1:Ntracks
    Frames{i}=linspace(1,1000,1000);
    Trace{i}(FrameArray{i}(1):FrameArray{i}(end),1)=zeros(size(TrackCell{i},1),1);
    Trace{i}(FrameArray{i}(1):FrameArray{i}(end),1)=TrackCell{i}(:,11);
    Trace{i}(isnan(Trace{i}))=0;
    Inz=find(Trace{i}>0);
    Frames{i}=Frames{i}(Inz);
    Trace{i}=Trace{i}(Inz);
end

%% plot
TraceVal=[];

for i=1:Ntracks
    hold on
    plot(Frames{i},Trace{i},'LineWidth',2);
    TraceVal=[TraceVal;Trace{i}];
end
    hold off
    title('Integrated Intensity vs. Frame number for Tus Spots');
    ylabel('Integrated Intensity (-)');
    xlabel('Frame (-)');
    set(gca,'FontSize',16);
    axis([0 200 0.01 80]);
    
%% plot of single trace
% ExTrace=nonzeros(Trace{1});
% Lt1=length(ExTrace);
% ExTrace=[ExTrace;Bkg.mean(Lt1:200)];
% plot(ExTrace)
% axis([0 200 200 400])
% xlabel('Slice (-)');
% ylabel('Mean Intensity (-)');
% title('Single Trace Intensity vs. Frames')
% set(gca,'FontSize',16);

%% save
SMCResult.utrack=utrack;
SMCResult.TraceVal=TraceVal;
SMCResult.TrackCell=TrackCell;

save(strcat(init.datapath,'/SMCResult',num2str(TraceNumber)),'SMCResult');





    