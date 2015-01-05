function [DivClicks,Division]=Division_GetCyclePoints(initval,DivClicks,Division,kymoprops, kymo_BF)
%Detect (manually) start and end points of a bacterium;
%that may consist of one or two spots but is considered to represent one
%replication cycle

%click cycle
stopit=0;
c=1;
clr='r'


while stopit==0;
%subplot(1,2,2); hold off
c=c+1;

%1) %pick living spots
%Find and analyse the living bacteria on this time
buf = struct2cell(DivClicks);
fate=squeeze(buf(2,:,:));
[~,Nbac]=size(DivClicks);  %number of bacteria
pcl=find(strcmp(fate, 'alive'));  %find indices of present bacteria  
PresentDivs=DivClicks(pcl) ;                    %...and their names

currentxposses=[]; currentyposses=[];    %for plotting purposes

if length(pcl)==0, 
    stopit=1;
else
  figure;
end

GoTonextGen=0;
for idx=1:length(pcl)
    bacno=pcl(idx); % 
   %idx: index of DivClicks in 'alive' DivClicks cycles
   %bacno: original index of DivClicks in database
ThisBac=PresentDivs(idx)

%show spot zoom right, total left

x1_abs=[ThisBac.PosClick.firstleft ThisBac.PosClick.firstright]; 
y1_abs=[ThisBac.PosClick.firstfr ThisBac.PosClick.firstfr];

%values used for viewing pane
zoomstart=ceil(ThisBac.PosClick.firstfr);
maxfr=min([kymoprops.duration ceil(zoomstart+kymoprops.zoom)]);
view_hor=initval.kymolength;
view_vert=maxfr-zoomstart;

%pcolor(kymo_BF(ceil(zoomstart):maxfr,:)); shading flat; colormap hot; title(titl); hold on;
x1_rel=x1_abs;      %coordinates relative to plotted kymograph section
y1_rel=[1 1];
[~,kl]=size(kymo_BF);

if GoTonextGen  %This ensures that all remaining 'live' bacteria will be skipped!
x2_rel=[-1 -1 -1];
y2_rel=[2 2 2];
but=[1 1 1];
else
titl=strcat('Zoomfrom',num2str(zoomstart),'to',num2str(maxfr));
P_Color(kymo_BF(zoomstart:maxfr,:),view_hor,view_vert,'hot'); title(titl); hold on;
plot(x1_rel,y1_rel,'-o','MarkerFaceColor', clr,'MarkerSize', 10); hold on; pause(0.5);
cor=1;
while cor==1
%---------------------------------------------------------------------
%explicit loop to show points
x2_LR=[];
y2_LR=[];
but=[];
for i=1:2  %click LEFT and RIGHT
[x,y,bt] = ginput(1);
plot(x,y,'-o','MarkerFaceColor', clr,'MarkerSize', 10); hold on;
x2_LR=[x2_LR x];
y2_LR=[y2_LR y];
but=[but bt ];
end
x2_rel=[x2_LR(1) (x2_LR(1)+x2_LR(2))/2 x2_LR(2)]; %LEFT MID RIGHT
y2_rel=[y2_LR(1) (y2_LR(1)+y2_LR(2))/2 y2_LR(2)]; %LEFT MID RIGHT
if max(but)==3, cor=1;else cor =0; end
end
end

Xav=mean(x2_rel);  %these are used for messaging ; see below
Yav=mean(y2_rel);

[a,i1]=min(x2_rel);
[c,i3]=max(x2_rel);
i2=find(x2_rel~=a & x2_rel~=c); b=x2_rel(i2);

[x2_rel,ix]=sort(x2_rel);
y2_rel=y2_rel(ix);

x2_rel=round(x2_rel);
y2_rel=round(y2_rel);

x2_abs=x2_rel;
y2_abs=y2_rel+zoomstart-1;

xplot=[x1_rel(1) x2_rel x1_rel(2)];
yplot=[y1_rel(1) y2_rel y1_rel(2)];   %inverted 'u' shape for plotting
if ~GoTonextGen & Xav<kymoprops.width
hold on; plot(xplot,yplot,'-ro','MarkerFaceColor', 'w');
end
pause(0.3);
close(gcf);

%1) update present cluster props %%%%%%%%%%%%%%%%%%%%%%%%%%%
xleft =min([x2_abs(1) kymoprops.width]);  %cut at end
xmid  =min([x2_abs(2) kymoprops.width]);  %cut at end
xright=min([x2_abs(3) kymoprops.width]);  %cut at end

yleft =min([y2_abs(1) kymoprops.duration]); %cut at top
ymid  =min([y2_abs(2) kymoprops.duration]); %cut at top
yright=min([y2_abs(3) kymoprops.duration]); %cut at top

DivClicks(bacno).PosClick.lastfr=ymid;
DivClicks(bacno).PosClick.lastleft=xleft;
DivClicks(bacno).PosClick.lastright=xright;

%Label current cycle; based on clicking %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Standard continuation-----------------

DivClicks(bacno).fate='divided';
GoTonextGen=0;
%UNLESS:
%1) %exit, go to next 'alive'   
if     Xav>kymoprops.width 
DivClicks(bacno).fate='exit'; 
GoTonextGen=0;
disp('exit')
end
%2) %exit; skip remaining set of 'alive' bacteria
if  Xav<0 
DivClicks(bacno).fate='exit';
GoTonextGen=1;
end

DivClicks(bacno).fate
%2) If division,Create two new clusters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp (DivClicks(bacno).fate,'divided')
if xmid<kymoprops.width & ymid<kymoprops.duration
%make a 'firstborn' offspring
DivClicks(Nbac+1).name=double(2*DivClicks(bacno).name);
DivClicks(Nbac+1).fate='alive';
DivClicks(Nbac+1).PosClick.firstfr=ymid;
DivClicks(Nbac+1).PosClick.firstleft=xleft;
DivClicks(Nbac+1).PosClick.firstright=xmid;
Division(Nbac+1).PosKyTrac.frames=ymid;
Division(Nbac+1).PosKyTrac.left=xleft;
Division(Nbac+1).PosKyTrac.right=xmid;
Nbac=Nbac+1;
end
if xright<kymoprops.width & ymid<kymoprops.duration
%make a 'secondborn(right)' offspring
DivClicks(Nbac+1).name=double(2*DivClicks(bacno).name+1);
DivClicks(Nbac+1).fate='alive';
Division(Nbac+1).PosKyTrac.frames=ymid;
Division(Nbac+1).PosKyTrac.left=xmid;
Division(Nbac+1).PosKyTrac.right=xright;
DivClicks(Nbac+1).PosClick.firstfr=ymid;
DivClicks(Nbac+1).PosClick.firstleft=xmid;
DivClicks(Nbac+1).PosClick.firstright=xright;
Nbac=Nbac+1;
end
end

%%%plot menu 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%absolute positions for plotting purposes (left-mid-right)
currentxposses=[currentxposses ; x1_abs(1) x2_abs x1_abs(2)];
currentyposses=[currentyposses ; y1_abs(1) y2_abs y1_abs(2)];


titl=strcat('Zoomfrom',num2str(zoomstart),'to',num2str(zoomstart+view_vert));
hold off ;  
%pcolor(kymo_BF(ceil(zoomstart):maxfr,:)); shading flat; colormap hot; title(titl); hold on;
P_Color(kymo_BF(zoomstart:maxfr,:),view_hor,view_vert,'hot'); hold on;;  colormap hot ;title(titl); hold on;
for p=1:i
%xplot=currentxposses(:,p);  yplot=currentyposses(:,p)-zoomstart; 
if min(xplot)<kymoprops.width & min(yplot)<kymoprops.duration
plot(xplot,yplot,'-o','MarkerFaceColor', clr);
end
end
hold on;
end
close(gcf);
outname=strcat(initval.basepath,kymoprops.WorkspaceOutName)
save(outname, 'Division', 'DivClicks','-append');
end

%Clean up
Bax = struct2cell(DivClicks);
fate=squeeze(Bax(2,:,:));   
subset=find(strcmp(fate, 'nonexistent')~=1);  %keep indices of present bacteria
DivClicks=DivClicks(subset);
Division=Division(subset);
DivClicks
