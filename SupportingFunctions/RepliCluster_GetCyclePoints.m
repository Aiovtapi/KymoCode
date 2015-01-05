function [ReplicationCluster,RepClicks]=RepliCluster_GetCyclePoints(chs,ReplicationCluster,RepClicks, kymoprops, kymo_FL,initval,actions)
%Detect (manually) start and end points of a replication spot 'cluster';
%that may consist of one or two spots but is considered to represent one
%replication cycle

 
 stopit=0;
 c=1;
 close all;


while stopit==0;   
%subplot(1,2,2); hold off
c=c+1;

%1) %pick living spots

%Find and analyse the living bacteria on this time
buf = struct2cell(RepClicks);
[~,Nrep]=size(RepClicks);  %current cycles stored
fate=squeeze(buf(2,:,:));   
pcl=find(strcmp(fate, 'present'));  %find indices of present bacteria  
PresentClusters=RepClicks(pcl) ;                    %...and their names
%number of 'alive' bacteria
currentxposses=[]; currentyposses=[];    %for plotting purposes
if length(pcl)==0, 
stopit=1;
else
figure;
end

GoTonextGen=0;


%Loop trough all the present bacteria -------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(pcl)
clusterno=pcl(i);

%show spot zoom right, total left
zoomstart=PresentClusters(1).PosClick.firstframe;
pos=PresentClusters(i).PosClick.firstpos;
%values used for viewing pane
maxfr=min([kymoprops.duration ceil(zoomstart+kymoprops.zoom)]);
view_hor=kymoprops.width;
view_vert=maxfr-zoomstart;

%plot the present zoomed part of the kymograph-------------------
%include overlay of division points
titl=strcat('Zoomfrom',num2str(zoomstart),'to',num2str(maxfr));
zoomkymo=kymo_FL(zoomstart:maxfr,:);
if actions.eliminate_growth, [zoomkymo,oripos]=Processing_Straighten_Growth(zoomkymo,initval);  
end
% close(gcf);
scrsz = get(0,'ScreenSize');
%figure('Position',[20 100 scrsz(3)/1.1 scrsz(4)/1.5]);
P_Color(zoomkymo,view_hor,view_vert,'jet'); title(titl); hold on;
plot(pos, 1,'o','MarkerFaceColor', 'k', 'MarkerSize', 12, 'MarkerEdgeColor','r'); hold on;

[lc,~]=size(currentxposses);
if lc>1;
    for p=1:lc
    xiplot=currentxposses(p,:);  yiplot=currentyposses(p,:)-zoomstart; 
    plot([xiplot(1) xiplot(3)],[yiplot(1) yiplot(3)],'w-o','MarkerFaceColor', 'w'); hold on;
    plot(xiplot(2:4),yiplot(2:4),'w-o','MarkerFaceColor', 'w'); hold on;
    %pause(0.01);
    end
end
pause(0.01);

%---------------------------------------------------------------
if GoTonextGen
xi=[1 1 1]'*1.2*kymoprops.width;
yi=[1 1 1]'*1.2*kymoprops.duration;
but=[1 1 1]';

else
cor=1;
while cor==1
    xi=[];  yi=[];
    for j=1:3 %click three main points; correct right-click
    [x,y,but] = ginput(1); 
    xi=[xi x]; yi=[yi y];
    plot(xi,yi, 'wo'); hold on;
    end
if max(but)==3, cor=1;else cor =0;end
if yi(1)> kymoprops.zoom, GoTonextGen=1, end  %skip remainders
end
end

close all


%figure;
%%%%%%sort points: left x early point  right x
[B,idx]=sort(yi); xi=xi(idx); yi=yi(idx);  %lowest time first
xi2=xi(2:3); yi2=yi(2:3); ;[B,idx]=sort(xi2); xi2=xi2(idx); yi2=yi2(idx);
xi=[pos xi2(1) xi(1) xi2(2) ]; 
yi=round([1 yi2(1) yi(1) yi2(2)]);

%%%absolute positions for plotting purposes
currentxposses=[currentxposses ; xi];
currentyposses=[currentyposses ; yi+zoomstart-1];


%%%plot menu 2 %%%%%%%%%%%%5
titl=strcat('Zoomfrom',num2str(zoomstart),'to',num2str(maxfr));
zoomkymo=kymo_FL(zoomstart:maxfr,:);

if actions.eliminate_growth, [zoomkymo,oripos]=Processing_Straighten_Growth(zoomkymo,initval);  
end
%close(gcf);
% scrsz = get(0,'ScreenSize');
% figure('Position',[20 100 scrsz(3)/1.1 scrsz(4)/1.5]);
P_Color(zoomkymo,view_hor,view_vert,'jet'); title(titl); hold on;
for p=1:i
xiplot=currentxposses(p,:);  yiplot=currentyposses(p,:)-zoomstart; 
plot([xiplot(1) xiplot(3)],[yiplot(1) yiplot(3)],'wo-','MarkerFaceColor', 'w','MarkerSize', 8,'LineWidth',2); hold on;
plot(xiplot(2:4),yiplot(2:4),'wo-','MarkerFaceColor', 'w','MarkerSize', 8,'LineWidth',2); 
hold on;
pause(0.01);
end

hold on;
%%%%%%%%%%%%%%%%%%%%%%%5


if actions.eliminate_growth  %coorect for the image deformation
    xi=xi.*2.^((yiplot-1)/initval.estimateddoublingtime);
end

%1) update present cluster props %%%%%%%%%%%%%%%%%%%%%%%%%%%
xstart=PresentClusters(i).PosClick.firstpos;    %current start pos
ystart=PresentClusters(i).PosClick.firstframe; %current start time
yend=min([yi(3)+zoomstart kymoprops.duration]); %cut at top
xend=min([xi(3) kymoprops.width]);  %cut at end

%2) Label end of current cluster--------------------------------------------
RepClicks(clusterno).fate='disassembled';
if     xend==kymoprops.width
RepClicks(clusterno).fate='exit';
end
if     xend<0
    RepClicks(clusterno).fate='exit';
    GoTonextGen=1;
end

%--------------------------------------------

%3) Set last positions current cluster----------------------
RepClicks(clusterno).PosClick.lastframe=round(yend);
RepClicks(clusterno).PosClick.lastpos=round(xend);

%4) make a first approximation of 'tracked' positions----------
life=yend-ystart+1;
bf=[ystart yend]; ystart=min(bf); yend=max(bf); %just to be sure
poss=linspace(xstart,xend,life);
frs=[ystart:1:yend];
ReplicationCluster(clusterno).PosKyTracCom.clickpos=poss;
ReplicationCluster(clusterno).PosKyTracCom.frames=frs;


if ~GoTonextGen
    %2) Create two new clusters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yleft=min([yi(2)+zoomstart kymoprops.duration]); %cut at top
    xleft=min([xi(2) kymoprops.width]);  %cut at end
    yright=min([yi(4)+zoomstart kymoprops.duration]); %cut at top
    xright=min([xi(4) kymoprops.width]);  %cut at end
    if xleft<kymoprops.width & yleft<kymoprops.duration
    %make a 'firstborn' offspring
    RepClicks(Nrep+1).name=double(2*RepClicks(clusterno).name);
    RepClicks(Nrep+1).PosClick.firstframe=yleft;
    RepClicks(Nrep+1).PosClick.firstpos=xleft;
    RepClicks(Nrep+1).fate='present';
    Nrep=Nrep+1;
    end
    if xright<kymoprops.width & yright<kymoprops.duration
    %make a 'secondborn(right)' offspring
    RepClicks(Nrep+1).name=double(2*RepClicks(clusterno).name+1);
    RepClicks(Nrep+1).PosClick.firstframe=yright;
    RepClicks(Nrep+1).PosClick.firstpos=xright;
    RepClicks(Nrep+1).fate='present';
    Nrep=Nrep+1;
end
end
end
h=gcf;
%print(h, '-dpng', '-r600',strcat(initval.basepath,initval.FiguresFolder,'ManualReplicationClicking/Channel',int2str(chs),'/InitializationStartPositions_Ch',int2str(chs)),'_Zoomfrom',num2str(zoomstart),'To',num2str(maxfr));
print(h, '-dpng', '-r150',strcat(initval.basepath,initval.FiguresFolder,'ManualReplicationClicking/Channel',int2str(chs),'/ReplicationClickingPositions_Ch',int2str(chs),'_ZoomFrom_',num2str(zoomstart),'To',num2str(maxfr) ))

% h=gcf;
% print(h, '-dpng', '-r600','/Volumes/LittleMonster/PhD_Docs/20130828_DnaN_Experiments/DnaN_Tiffs_Series1002 37C/Figures/RepCycleClickingTest')
end
%Clean up
Bax = struct2cell(RepClicks);
fate=squeeze(Bax(2,:,:));   
subset=find(strcmp(fate, 'nonexistent')~=1);  %keep indices of present bacteria
RepClicks=RepClicks(subset);
ReplicationCluster=ReplicationCluster(subset);
end