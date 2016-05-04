function RepliCluster00_TrackandcleanQuick(exp,ColourIdx)
%ReplicationCluster
%-------------------------------------------------------------------------
if nargin<1, exp='001_DnaN_TUS_dif_30122014_TUSsignal';end 

initval=A001_Images_Set_Experiment(exp);
chans=initval.channelno;
for ch=1:chans
close all
display(strcat('chans to go = ',num2str(chans-ch)));
Channelpath=char(strcat(initval.basepath,initval.nms{ch}(ColourIdx),'.mat'));
load(Channelpath);

WorkspaceOutName=char(initval.nms{ch}(ColourIdx)); % ugly way to update entries
    
[r,c]=size(kymo_FL);
kymoprops.width=c;
kymoprops.duration=r;
kymoprops.zoom=70;  %used for clicking

actions.cleandatabase=1;
actions.dotracking=1;
actions.cleantracking=1;


if actions.cleandatabase==1;
    [ReplicationCluster,RepClicks]=RepliCluster_Cleaning(ReplicationCluster,RepClicks);
end 
if actions.dotracking==1
    ReplicationCluster=RepliCluster_DoTrackingQuick(ReplicationCluster,kymoprops,kymo_FL,initval);
end


%Plotting menu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,le]=size(ReplicationCluster);
close all;
pcolor(kymo_FL); shading flat; colormap hot; title('overlay'); hold on;
for i=1:le  
    repno=i;
    pos=ReplicationCluster(i).PosKyTracCom.trackpos; 
    ti=ReplicationCluster(i).PosKyTracCom.frames;
    %plot(pos+1,ti, 'b-'); hold on
    plot(pos+0.5,ti+0.5, 'bo-');
    title(strcat('Edges Analysis ',WorkspaceOutName));    
end
 pause(5);
 
FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'TrackAndCleanQuick/'));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'TrackAndCleanQuick/'));
end

h=gcf;
print(h, '-dpng', '-r150',strcat(initval.basepath,initval.FiguresFolder,'TrackAndCleanQuick/',WorkspaceOutName))

%  h=gcf;
%  print(h, '-dpdf', '-r600','Testing')

outname=strcat(initval.basepath,WorkspaceOutName);
save(outname, 'initval', 'RepClicks', 'ReplicationCluster',  '-append');
end

