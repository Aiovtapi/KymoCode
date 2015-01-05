function dum=Binfiller(Xdata,Ydata,binax,initval,extras);
% This function allocates data to a series of bins and averages the result
% JacobKers 2013
dum=1;
scatwidth=0.7;
Lbins=length(binax);
minbin=binax(1);
maxbin=binax(end);
ld=length(Ydata);
maxno=ld/Lbins;
binsize=(maxbin-minbin)/(Lbins-1);

bincollector=zeros(Lbins,ld)*NaN;  %define bins
bincounts=zeros(Lbins,1); %keep track how many datapoints are added per bin

binidx=ceil((Xdata-minbin)/(maxbin-minbin)*Lbins); %bin number for each data  point
% los=find(binidx<minbin); binidx(los)=minbin;        %collect outliers in extrem bins
% his=find(binidx>maxbin); binidx(his)=maxbin;

los=find(binidx<1); binidx(los)=1;        %collect outliers in extrem bins
his=find(binidx>Lbins); binidx(his)=Lbins;
 

%fill bins--------------------------------------------------------------
for i=1:ld   
ix=binidx(i);
bincounts(ix)=bincounts(ix)+1;  %adjust bin counter
bincollector(ix,bincounts(ix))=Ydata(i);  %put point of curve in proper bin
end

%process collected bin data-------------
filledbinsidx=find(bincounts~=0);  %keep only filled bins
bincounts=bincounts(filledbinsidx);
binaxnw=binax(filledbinsidx);
Lbins2=length(binaxnw);
bincollector=bincollector(filledbinsidx,:);

%Define a 'scatter' axis and an axis for each individual entry in the bins
scatax2D= scatwidth*binsize*(rand(Lbins2,ld)-0.5);
binax2D=repmat(binaxnw',1,ld);

%clean the peak data from outliers and get averages-----------------
av=zeros(Lbins2,1);
sig=zeros(Lbins2,1);

for i=1:Lbins2    
sel=~isnan(bincollector(i,:));
onebindata=bincollector(i,sel);
%Detect 'inliers' and refine selection--------------------
[flag,cleandata]=Outlier_Flag(onebindata,3,0.8,'positive',0,20);
%data,tolerance,sigchange,how,sho, binz
sel2=find(flag==1);
%Get averages and standard deviation of final selection
av(i)=nanmean(onebindata(sel2));
sig(i)=nanstd(onebindata(sel2));
end
%Finally, upper and lower bounds (one sigma)-------------------------------
lo=av-1*sig;
hi=av+1*sig;

%Plot menu--------------------------------------------------------------
%1a)  total&peaks per bacterium per time point-----------------------------
h=figure;
h1 = plot(binax2D+scatax2D,bincollector, 'o', 'color' ,[0.4 0.4 0.4], 'MarkerSize', 2); 
h1=h1(1);
hold on;
%plot(binax2D+scatax2D,bincollector, 'bo', 'MarkerSize', 1.5); hold on;
title(strcat('Decrease in fluorescence as function of time:',extras.remark),'fontsize', 16, 'fontweight', 'bold');
xlabel('Frames (-)', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Intensity per cell length (-)', 'fontsize', 16, 'fontweight', 'bold');
axis([minbin max(binaxnw) 0 max(Ydata)]);
set(gca, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%1b totals peaks av-lo-hi
h2 = plot(binaxnw,av,'--ro','LineWidth',3,...
                 'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on;
h3 = plot(binaxnw,lo,'b-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerSize',4); hold on;
plot(binaxnw,hi,'b-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerSize',4); hold on;
legend([h1 h2 h3],{'Binned data points', 'Average', 'Standard deviation'});
 

FileName_bleachingCurve = 'BleachingAfterCommencementExperiment';
lbl=strcat(initval.basepath,FileName_bleachingCurve);
print(h, '-dpdf', '-r600',lbl)           
hold on;           
%--------------------------------------------------------------------------

%Saving menu---------------------------------------------------------------
if 1
     %2) binned averages
     avdata=[binaxnw' lo av hi];
    filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',initval.expno,extras.remark,'_BinTotals.txt');
     dlmwrite(filname,avdata, 'delimiter', '\t');
     %2) individual points per bin
     x=(binax2D+scatax2D);
     x1=binax2D;
     binpoints=double([x1(:) x(:) bincollector(:)]);
    filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',initval.expno,extras.remark,'_BinSinglePoints.txt');
    format long;
     dlmwrite(filname,binpoints, 'delimiter', '\t','precision','%10.5f');
end

