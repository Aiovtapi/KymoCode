function B0011_CytoplasmicConcentration
%Load database, make binned curves
%JacobKers 2012----------------------------------
close all; clear all

%exp='TEST';
%exp='CM_DnaN-Dif-Gamma-ve-Position1_Series1';
exp='001_DnaN_TUS_dif_30122014_DnaNsignal';

actions.loaddatabase=1;         %default=1 (new analysis)

remark='_HFR_Exp1';  %This allows you to add comments to the filenames
titlelabel='First Higher Frame Rate experiment';
minsperframe=2.5;

%The name of the figure files that are printed to PDFs
FileName_Percentage='FluorescenceCell_Foci'; %used for drift correction
FileName_RepDivHistograms='DivRepTimes_Hist_Combined_37degrees';
FileName_FociCyto = 'Total_and_Foci_Increase_Initiation_Termination';

%Settings for building histograms
maxbins=50;
scatwidth=0.8;
countsperlabel=3600;
maxax=1E5/countsperlabel;  %Sets scale axis graphs
timeax='relative';      %for cycle
maxbintime=1.2;         %for div/rep histogram
minbintime=-0.2;
absmaxbintime=200;      %for div/rep histogram

initval=A001_Images_Set_Experiment(exp);

%load the databases, build arrays--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
[~,chan_no]=size(S);
[good, bad]=Processing_Measure_GoodDatabase(S)
maxbac= good.birthdivedgesuser;
bincollector=zeros(maxbins,2*maxbac)*NaN;  %array for collecting averaged data
%-------------------------------------

%we assume 100x2.5 mins will fit all bacteria and each bacterium will add
%two datapoints per bin
bincounts=zeros(maxbins,1); %keep track how many datapoints are added per bin
binax=minbintime+(maxbintime-minbintime)/maxbins*[0.5:maxbins-0.5];  %time axis 
absbinax=absmaxbintime/maxbins*[0.5:maxbins-0.5];
totalbaccount=0;
alldivtime=[];
allreptime=[];
AllBax_times=[];
AllBax_totals=[];

%%
for ch=1:chan_no  %for each channel
chan_no-ch;
Div=S(ch).channels.AutoDivision;
Rep=S(ch).channels.ReplicationCluster;

kymo_FL=S(ch).channels.kymo_FL;
stripmov_FL=S(ch).channels.chanstk_FL;
[~,bacno]=size(Div);
for j=1:bacno  %for each bacterium 
    % bacno-bc
    bacname=S(ch).channels.RepClicks(j).name;
    ThisRep=Rep(j);         %Current replication
    ThisBac=Div(j);
    frs= ThisBac.edges.frs; 
    %check number of conditions   
    Okay=Check_Bacterial_Quality(ThisBac,ch,j,S);
     if Okay       
        totalbaccount  =totalbaccount+1;  
         %Select the fluorescent values-------------------------------------
        %0) Select the time axis (normalized or not)
        switch timeax
            case 'absolute', %time in minutes, starts at 0 per bacterium
                frames=minsperframe*S(ch).channels.ReplicationCluster(j).PosKyTracCom.frames_ext;
                frames=frames-frames(1);  
            case 'relative', %time in scaled units, 0 at init, 1 at ter
                frames=S(ch).channels.ReplicationCluster(j).Cycle.normtime_ext';
        end
        %------------------------------
        %1) Select time axis, data of interest 
        abs_frames=S(ch).channels.ReplicationCluster(j). FluoPropsGen.frs_ext'; 
        fluocounts=S(ch).channels.ReplicationCluster(j).FluoPropsGen.fluosumspinelevel_ext;

        %In this section, data is corrected for the expected bleaching ------
        labels=Convert_fluorescence_to_labelcounts(fluocounts,abs_frames);       
        dimers=Convert_labelcounts_to_dimers(labels);  
        nMol=Convert_dimers_to_concentration(dimers);  
        %-------------------------------------------------------------------------
        %In this section, we concatenate bacterium data for later processing
        lp=length(frames);
        idxes=[1:lp]';
        AllBax_times=[AllBax_times ; [idxes frames]];
        AllBax_totals=[AllBax_totals ;[idxes nMol']];
        end
    end
end


%%
%-----------------------------------------------------------------------
%Run Binfillers: In this section, data is redivided over bins. 
%We choose to have one data point per bacterium per bin; Therefore the data is resampled.
binresults=Binfiller_General_IntPol(AllBax_times,AllBax_totals,binax,0);
%-----------------------------------------------------------------------


%---------------------------------------------------
%In this section, we make  stacked histograms of the data points over a period of
%normalized time, for example from t=0.2 to 0.8, where t=0 indicates
%initiation and t=1 indicates termination. The aim is to see if outliers
%affect our average & sigma too much. Note : this action only makes sense
%if we use relative times hence the condition


if strcmp(timeax,'relative');      %for cycle
timeslots=[-0.1:0.2:1.1];  %these are begin and end points of the relative time intervals
nMol_bins=linspace(0,200,20);
dum=Proc_Make_histograms_binresults(binresults,timeslots, nMol_bins);
end
%-------------------------------------------------------------






%General plot and save menu--------------------------------------------------------
%1a) Individual points and averages per bin--------------------------------------------------------------
%for plotting use, extract data-----------------
ax2D=binresults.binaxis_all;
axscat2D=binresults.binaxis_all_scat;
all_cell=binresults.binvalues_all;

loall=binresults.losig;
avall=binresults.av;
hiall=binresults.hisig;
%----------------------------


% %%------------------------HERE WE PLOT THE NUMBER OF MOLECULES (REL/ABS) IN CYTOPLASM
% %%------------------------AS WELL AS INSIDE THE FOCI
h=figure;
h1 = plot(axscat2D,all_cell, 'bo', 'MarkerSize', 1.5); hold on;
h1=h1(1);
title(strcat('Cytoplasmic Concentration , N=',num2str(totalbaccount)), 'fontsize', 16, 'fontweight', 'bold')
xlabel('Norm. Time Init-Ter (-)', 'fontsize', 16, 'fontweight', 'bold');
ylabel('concentration, nM', 'fontsize', 16, 'fontweight', 'bold');
%ylabel('Number of \beta_2 clamps (-)', 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%axis([minbintime max(binax) 0 maxax]);
%1b totals peaks av-lo-hi

h3 = plot(binax,avall,'--bo','LineWidth',2,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on;
h4 = plot(binax,loall,'b-','LineWidth',3,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,hiall,'b-','LineWidth',3,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
            
hleg1 = legend([h1 h3],{'Binned data points in cell','Average intensity in cell'});
set(hleg1,'Location','NorthWest')
      
lbl=strcat(initval.basepath,FileName_FociCyto);
%print(h, '-dpdf', '-r600',lbl)
%print(h, '-dpdf', '-r600','FociCyto_DnaN')

% dlmwrite('CytoplasmicConcentration_WT_AVG.txt', [binax' avall], 'delimiter', '\t','precision', 3)
% dlmwrite('CytoplasmicConcentration_WT_LOW.txt', [binax'  loall], 'delimiter', '\t','precision', 3)
% dlmwrite('CytoplasmicConcentration_WT_HIGH.txt', [binax' hiall], 'delimiter', '\t','precision', 3)


% dlmwrite('CytoplasmicConcentration_GammaMinus_AVG.txt', [binax' avall], 'delimiter', '\t','precision', 3)
% dlmwrite('CytoplasmicConcentration_GammaMinus_LOW.txt', [binax'  loall], 'delimiter', '\t','precision', 3)
% dlmwrite('CytoplasmicConcentration_GammaMinus_HIGH.txt', [binax' hiall], 'delimiter', '\t','precision', 3)
%    %dlmwrite('Hist_Data_WT.txt', [binCtrs' prob'], 'delimiter', '\t','precision', 3)


%Saving menu---------------------------------------------------------------
     %2) binned averages
     avdata=[binax'-minsperframe/2 loall avall hiall];
    filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',exp,remark,'_SpotsAndTotals.txt');
     dlmwrite(filname,avdata, 'delimiter', '\t');
     
     %2) individual points per bin
     binpoints=[ax2D-minsperframe/2 axscat2D-minsperframe/2 all_cell];
    filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',exp,remark,'_BinPoints.txt');
     dlmwrite(filname,binpoints, 'delimiter', '\t');
     
     
  function Okay=Check_Bacterial_Quality(ThisBac,ch,j,S);
     ok1=strcmp(ThisBac.birthtype, 'OK');    %birth ok
     ok2= strcmp(ThisBac.divtype, 'OK');    ;%division ok
     ok3=ThisBac.edges.edgesok;  %edges ok
     frs= ThisBac.edges.frs; 
     ok4=(length(frs)>0); 
     ok5=strcmp(S(ch).channels.RepClicks(j).fate, 'disassembled');
     ok6=S(ch).channels.AutoDivision(j).accepted;  %manual accept/reject
     Okay=ok1&ok2&ok3&ok4&ok5&ok6;
 
     
     
     
function correctedlabels=Convert_fluorescence_to_labelcounts(fluo,frames);
%function corrects for bleaching; JacobKers 2013
countsperbetaclamp=3600;
%This is the conversion unit when no bleaching is present
decreasefactor=2.1;  
%This is a factor describing how much the fluorescence
%of bacteria will have dropped when bleaching and expression balannce in a 'steady-state' signal 
decay=30; 
%This is the time constant of bleaching in frame units, assuming
%exponential decay to the steady-state level
%Above factors must be determined elsewhere.
if nargin<2  %if TEST MODE
    close all;
    frames=[10:60];
    labels=0*frames+100;   %constant 100 labels
    bleachfactorperframe=((decreasefactor-1)*exp(-frames/decay)+1)/decreasefactor;
    bleachfactorperframe=((decreasefactor-1)*exp(-frames/decay)+1)/decreasefactor;
    fluo=countsperbetaclamp*labels.* bleachfactorperframe;     
end
%-------------------------------------------------------------------
correctfactor=1./(((decreasefactor-1)*exp(-frames/decay)+1)/decreasefactor)';
uncorrectedlabels=fluo/countsperbetaclamp;

correctedlabels=correctfactor.*uncorrectedlabels;
if nargin<2 %TEST MODE
     plot(frames,labels,'-'); hold on;
     plot(frames,uncorrectedlabels,'-o');
     plot(frames,correctedlabels,'ro');
     legend('original labels', 'calculated labels, no bleachcorrection', 'calculated labels, bleachcorrected' );
end

function dimers=Convert_labelcounts_to_dimers(labels);
    dimers=labels/2;
    
function muMol=Convert_dimers_to_concentration(dimers); 
      %This function returns an estimate of the cytoplasmic concentration
      %in micromolair
      %using the median sum of dimers along the perpendicular direction of a
      %bacterium, summed per unit pixel length
      p2m=160E-9;    %unitpixel size in meters
      br=0.35E-6;    %bacterium_radius in meters
      vol=p2m*pi*(br^2)*1000;  %volume in liters assuming circular radius
      muMol=dimers/6.02214129E23/vol*1E9;
 
      
      function dum=Proc_Make_histograms_binresults(binresults,slots, binss);
%Make histograms of binned results. JacobKers 2013
%---------------------------------------------------
%In this section, we make  stacked histograms of the data points over a period of
%normalized time, for example from t=0.2 to 0.8, where t=0 indicates
%initiation and t=1 indicates termination. The aim is to see if outliers
%affect our average & sigma too much. Note : this action only makes sense
%if we use relative times 
%example slots=[0.1  0.5  0.9];  %these are begin and end points of the relative time intervals
%Contents of 'results'
% binresults.losig;          %one sigma below the average, per bin
% binresults.av;             % the average, per bin
% binresults.hisig;          %one sigma above the average, per bin
% binresults.binaxis;   %values of the non-empty bins, per bin
% binresults.binaxis_all; %value of  bin, per point
% binresults.binaxis_all_scat; %random number to spread bin points (for plotting purposes)
% binresults.binvalues_all; %values themselves
% binresults.bincounts=traceno;

if nargin<2  %TEST modus with 50 fake cels
    close all; 
    slots=[0.1 0.3 0.6 0.9];
    binresults.binaxis_all=repmat(linspace(-0.1,  1.1,100),1,50);
    binresults.binvalues_all=repmat(3+randn(1,100),1,50);
    binz=40;
    lobin=0;
    %hibin=max(binresults.binvalues_all);
    hibin=200;
    binss=linspace(lobin,hibin,binz);
end


slotno=length(slots)-1;
for i=1:slotno
    sel=find(binresults.binaxis_all>slots(i) & binresults.binaxis_all<slots(i+1));
    buf=binresults.binvalues_all(sel);
    hist_slot=hist(buf,binss);
    subplot(slotno,1,i);
    bar(binss,hist_slot); hold on;
    text(binss(1)+0.6*(binss(end)-binss(1)),0.8*max(hist_slot),strcat('From',num2str(slots(i)),'to',num2str(slots(i+1))));
    if i==1, title('distribution per relative time {init-ter}'), end
    if i==slotno; xlabel('Label no., a.u.'),    end
    axis tight
end
dum=1;

function binresults=Binfiller_General_IntPol(Xdata2D,Ydata2D,binax,plotit);
% This function allocates data to a series of bins and averages the result;

% If the number of points in the X-axis is small (< bins) and normalized, direct filling is irregular. 
% In this version, the data is first interpolated to ensure each bin is
% filled once and only once per trace
% Use this version if you are  not interested in counts, but in averages


%input: X data, Y data: 2D array of index, values. Xdata en Y data is assumed to be
%concatenated; X data is assumed to increase in time per sub-trace
%output:
    % binresults.losig=lo;          %one sigma below the average, per bin
    % binresults.av=av;             % the average, per bin
    % binresults.hisig=hi;          %one sigma above the average, per bin
    % binresults.binaxis=binaxnw;   %values of the non-empty bins, per bin
    % binresults.binaxis_all=binax2D(:); %value of  bin, per point
    % binresults.binaxis_all_scat=scatax2D(:); %random number to spread bin points (for plotting purposes)
    % binresults.bincounts=bincounts;
    
% JacobKers 2013
%------------------------------------------------------
test=0;  %set to 1 for demo
scatwidth=0.7;

%---------------------------------------------------
if test
    close all
    plotit=1; %set to 1 for viewing results
    minbin=-0.1;
    maxbin=1.1;
    bins=100;
    binax=linspace(minbin,maxbin,bins);
    binsize=binax(2)-binax(1);
    %fake data
    [Xdata2D,Ydata2D]=CreateFakeData;
end
%----------------------------------------------------
bins=length(binax);
binsize=binax(2)-binax(1);
minbin=binax(1);
maxbin=binax(end);

Idxdata=Xdata2D(:,1);
Xdata=Xdata2D(:,2);
Ydata=Ydata2D(:,2);

%remove outliers looking at the whole dataset
[flag,~]=Outlier_Flag(Ydata,3,0.8,'all',0,50);
sel=find(flag==1);
Xdata=Xdata(sel);
Ydata=Ydata(sel);
Idxdata=Idxdata(sel);
% -------------------------------

%First, measure the number of separate datasets by looking at 'negative time
%jumps'; define bincollector matrix
dif=Idxdata(2:end)-Idxdata(1:end-1);
sel=(find(dif<0)); 
traceno=length(sel)+1;
tracestart=([0 sel' length(Idxdata)])';
bins_L=length(binax);
bincollector=zeros(bins_L,traceno)*NaN;  %define bins

%Next, transform each dataset and fill bins
for i=1:traceno
    x=Xdata(tracestart(i)+1:tracestart(i+1))';
    y=Ydata(tracestart(i)+1:tracestart(i+1))';
    [xq,ixq,~]=unique(x);   %remove double entries from earlier steps
    yq=y(ixq);
    xip=binax';
    yip=interp1(xq,yq,xip,'linear', NaN);
    bincollector(:,i)=yip;
    if 0; 
    plot(x,y,'o-'); hold on;
    plot(xip,yip,'ro'); hold off;    
    [~]=ginput(1);
    end 
    dum=1;
end

%Define a 'scatter' axis and an axis for each individual entry in the bins
scatax2D= scatwidth*binsize*(rand(bins_L,traceno)-0.5);
binax2D=repmat(binax',1,traceno);

%clean the peak data from outliers and get averages-----------------
av=zeros(bins_L,1);
sig=zeros(bins_L,1);

for i=1:bins_L    
sel=~isnan(bincollector(i,:));
onebindata=bincollector(i,sel);
    if length(onebindata)>0
    %Detect 'inliers' and refine selection--------------------
    [flag,cleandata]=Outlier_Flag(onebindata,2,0.8,'positive',0,10);
    %data,tolerance,sigchange,how,sho, binz
    sel2=find(flag==1);
    nonsel=find(flag==0);
    bincollector(i,sel(nonsel))=NaN;
    
    %Get averages and standard deviation of final selection
    av(i)=nanmean(onebindata(sel2));
    sig(i)=nanstd(onebindata(sel2));
    else
    av(i)=0;
    sig(i)=0;
    end
end
%Finally, upper and lower bounds (one sigma)-------------------------------
lo=av-1*sig;
hi=av+1*sig;

binresults.losig=lo;          %one sigma below the average, per bin
binresults.av=av;             % the average, per bin
binresults.hisig=hi;          %one sigma above the average, per bin
binresults.binaxis=binax;   %values of the non-empty bins, per bin
binresults.binaxis_all=binax2D(:); %value of  bin, per point
binresults.binaxis_all_scat=binax2D(:)+scatax2D(:); %random number to spread bin points (for plotting purposes)
binresults.binvalues_all=bincollector(:);
binresults.bincounts=traceno;


%Plot menu--------------------------------------------------------------
%1a)  total&peaks per bacterium per time point
if plotit
figure;
plot(binax2D+scatax2D,bincollector, 'bo', 'MarkerSize', 1.5); hold on;
title(strcat('Binned intensities-'));
xlabel('frames');
ylabel('Intensity per length');
axis([minbin max(binax) 0 max(Ydata)]);

%1b totals peaks av-lo-hi
plot(binax,av,'--bs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,lo,'bo-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,hi,'bo-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on; 
%--------------------------------------------------------------------------
end

function [Xdata,Ydata]=CreateFakeData;
%This function generates traces with length varying randomly around an average
%input: mininimum value X, maximum value X.
%output: Ntrace of these concatenated, X and Y
%It shows the '0' and '1' gaps artefacts
%JacobKers 2013

Ntrace=200;
Avpointspertrace=20;
lox=-0.3;
hix=1.8;
loy=20;
hiy=40;
Xdata=[]; Ydata=[];

for i=1:Ntrace
tracelength=ceil(Avpointspertrace*(1+0.5*(rand(1,1)-0.5))); %50% plus min variation
x=linspace(lox,hix,tracelength)';
y=linspace(loy,hiy,tracelength)'+0.5*(hiy-loy)*(rand(1,tracelength)-0.5)';
if 0; 
    plot(x,y,'o-');
    %[~]=ginput(1);
end
idx=(1:tracelength)';
Xdata=[Xdata ; [idx x]]; 
Ydata=[Ydata ; [idx y]];
end

