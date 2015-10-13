function RepliCluster00_TrackandcleanQuick(exp)
%ReplicationCluster
%-------------------------------------------------------------------------
if nargin<1, exp='001_DnaN_TUS_dif_30122014_TUSsignal';end 

initval_ori=A001_Images_Set_Experiment(exp);
chans=length(initval_ori.nms);
for ch=1:chans
close all
display('chans to go');
chans-ch
Channelpath=char(strcat(initval_ori.basepath,initval_ori.nms(ch),'.mat'));
load(Channelpath);

WorkspaceOutName=char(initval_ori.nms(ch)); % ugly way to update entries
initval=initval_ori; %to ensure overrulling settings
    
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

