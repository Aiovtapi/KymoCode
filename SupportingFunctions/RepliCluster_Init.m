function [ReplicationCluster,RepClicks]=RepliCluster_Init(chs,kymo_FL,kymoprops, initval,actions);             
 

zoomkymo=kymo_FL(1:initval.startzoom,:);
if actions.eliminate_growth, 
    [zoomkymo,oripos]=Processing_Straighten_Growth(zoomkymo,initval);  %TEST
end

%-------------------------------------------------------------------------
titl=strcat('Initialization: Click Start positions, then right-click channel end position');
scrsz = get(0,'ScreenSize');
figure('Position',[20 100 scrsz(3)/1.1 scrsz(4)/1.5]);
P_Color(zoomkymo,initval.kymolength,initval.startzoom,'hot'); hold on; colormap jet ;title(titl); hold on;
 
posvector=[];
tinit=[];
stop=0;
while stop==0
    [x,y,but]=ginput(1);
    %set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
    plot(x,y,'ro','MarkerSize', 10,'LineWidth',4); hold on;
    
%     
%     plot(driftvector(:,1),'-bo','LineWidth',4,...
%                 'MarkerEdgeColor','b',...
%                 'MarkerFaceColor','w',...
%                 'MarkerSize',6);
%set(gca, 'YTick', [0 - 0.02 0.03]);

% %set(gca,'TickLength',[0.02 0.02]);
% ylabel('Position (px)', 'fontsize', 26, 'fontweight', 'bold');
% xlabel('Frames (-)', 'fontsize', 26, 'fontweight', 'bold');

    
    
    posvector=[posvector x];
    tinit=[tinit y];
    if but==3, stop=1; end
end
FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'ManualReplicationClicking/Channel',int2str(chs)));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'ManualReplicationClicking/Channel',int2str(chs)));
end

h=gcf;
print(h, '-dpng', '-r150',strcat(initval.basepath,initval.FiguresFolder,'ManualReplicationClicking/Channel',int2str(chs),'/InitializationStartPositions_Ch',int2str(chs)))
if actions.eliminate_growth  %coorect for the image deformation
    posvector=posvector.*2^((y-1)/initval.estimateddoublingtime);
end

le1=length(posvector)-1;    %starting number, 
le2=2^(ceil(log2(le1))); %2^n to keep labelling consistent
padd=posvector(le1+1)*[1:le2-le1];
paddinit=1*[1:le2-le1];
posvector=ceil([posvector(1:le1) padd]);
tinit=ceil([tinit(1:le1) paddinit]);

%lp=length(posvector);

for i=1:le2  %only initialize entries necesary for first steps
    RepClicks(i).name=le2+i-1;
    RepClicks(i).fate='present';
   
    RepClicks(i).PosClick.firstframe=tinit(i);
    RepClicks(i).PosClick.lastframe=tinit(i)+1;
    RepClicks(i).PosClick.firstpos=posvector(i);
    RepClicks(i).PosClick.lastpos=posvector(i);
    
    ReplicationCluster(i).PosKyTracCom.frames=tinit(i);   
    ReplicationCluster(i).PosKyTracCom.clickpos=posvector(i);
end
for j=le1+1:le2
     RepClicks(j).fate='exit';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%