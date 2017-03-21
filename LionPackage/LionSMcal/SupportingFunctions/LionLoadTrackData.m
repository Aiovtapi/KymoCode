function [TrackCell,FrameArray,Bkg] = LionLoadTrackData(utrack,Bkg)
% This function translates the utrack output data to create a single cell
% with trajectory information in rows, instead of columns. This is the
% format:
%
% For each track:
% 
% create a cell containing matrices [X Y Z AMP Xstd Ystd Zstd Ampstd BkgMean Bkgstd]

lengthlowerbound=4;
umperpx=0.1589;

Ntracks=size(utrack.tracks.tracksFinal,1);
TrackCell=cell(Ntracks,1);
FrameArray=cell(Ntracks,1);

for i=1:Ntracks
    Noricolumns=size(utrack.tracks.tracksFinal(i).tracksCoordAmpCG,2);
    for j=1:Noricolumns/8;
    TrackCell{i}(j,:)=utrack.tracks.tracksFinal(i).tracksCoordAmpCG(1+8*(j-1):8*j);
    TrackCell{i}(j,4)=(TrackCell{i}(j,4))*2^16-1; % Amp Convert to Photons & integrated intensity
    TrackCell{i}(j,8)=TrackCell{i}(j,8)*(2^16-1); % STD Amp Convert to Photons 
    end
end

for i=1:Ntracks
    % add the background integrated intensity
    FrameArray{i}=linspace(utrack.tracks.tracksFinal(i).seqOfEvents(1,1),...
        utrack.tracks.tracksFinal(i).seqOfEvents(2,1),utrack.tracks.tracksFinal(i).seqOfEvents(2,1)...
        -utrack.tracks.tracksFinal(i).seqOfEvents(1,1)+1)';
    
        TrackCell{i}(:,9)=Bkg.mean(FrameArray{i});
        TrackCell{i}(:,10)=Bkg.std(FrameArray{i});
        
        % integrated intensity (normalized).
        TrackCell{i}(:,11)=(TrackCell{i}(:,4))*TrackCell{i}(j,5)*TrackCell{i}(j,6)*2*pi/umperpx^2;
        
    if size(TrackCell{i},1)<lengthlowerbound
        TrackCell{i}=[];
        FrameArray{i}=[];
    end
end

TrackCell=TrackCell(~cellfun('isempty',TrackCell));
FrameArray=FrameArray(~cellfun('isempty',FrameArray));


end