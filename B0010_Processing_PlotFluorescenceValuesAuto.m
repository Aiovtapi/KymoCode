%Processing_PlotFluorescenceValuesAuto
%Load database, make binned curves
%JacobKers 2012----------------------------------
close all
clear all
%select the experiment database-----------------------------

exp='001_DnaN_TUS_dif_30122014_DnaNsignal_P';


actions.savedata=1;             %default=1 (new analysis)
actions.loaddatabase=1;         %default=1 (new analysis)
variants={'_all'};  %'_mother' 'M1D1' see line 96, allows different selections

remark='__HFR_Exp1';  %This allows you to add comments to the filenames
titlelabel='First Higher Frame Rate experiment';

minsperframe=5;
maxbins=50;


%The name of the figure files that are printed to PDFs
FileName_Percentage='FluorescenceCell_Foci'; %used for drift correction
FileName_RepDivHistograms='DivRepTimes_Hist_Combined_37degrees';
FileName_FociCyto = 'Total_and_Foci_Increase_Initiation_Termination';

scatwidth=0.8;
%countsperlabel=1;
countsperlabel=3600;
maxax=1E5/countsperlabel;  %Sets scale axis graphs
spotsanalysistype='1DGauss'; %'medianexcess' %'area summing';

%spotsanalysistype='1DGauss' %'area summing';
%spotsanalysistype='area summing' %'area summing';
%spotsanalysistype='medianexcess_suppressotherclusters' %'area summing';

switch spotsanalysistype
    case '1DGauss', sptyp=1;
    case 'peakvalestimate', sptyp=2;
    case 'medianexcess', sptyp=3;
    case 'medianexcess_suppressotherclusters', sptyp=4;
    case 'area summing', sptyp=5;
    case 'DnaX', sptyp=6;
end

switch 2
    case 1
    timeax='absolute'; 
    maxbintime=160;         %for cycle
    minbintime=0;
    absmaxbintime=160;      %for div/rep histogram
    case 2
    timeax='relative';      %for cycle
    maxbintime=1.2;         %for div/rep histogram
    minbintime=-0.2;
    absmaxbintime=200;      %for div/rep histogram
end



variant=char(variants);

initval=A001_Images_Set_Experiment(exp);
outname=strcat(initval.basepath,initval.outname);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------

[~,chan_no]=size(S);
[good, bad]=Processing_Measure_GoodDatabase(S)

%arrays for collecting averaged data
%-------------------------------------
maxbac=300;

bincollector_all_labels=zeros(maxbins,2*maxbac)*NaN;  
bincollector_peak_labels=zeros(maxbins,2*maxbac)*NaN;
bincollector_perc_labels=zeros(maxbins,2*maxbac)*NaN;
%we assume 100x2.5 mins will fit all bacteria and each bacterium will add
%two datapoints per bin
bincounts=zeros(maxbins,1); %keep track how many datapoints are added per bin
binax=minbintime+(maxbintime-minbintime)/maxbins*[0.5:maxbins-0.5];  %time axis 
absbinax=absmaxbintime/maxbins*[0.5:maxbins-0.5];
%column 1: average of freelabels
%column two: standard deviation
totalbaccount=0;
alldivtime=[];
allreptime=[];
%-------------- -----------------------
tic
AllBax_times=[];
AllBax_totals=[];
AllBax_spots=[];
AllBax_perc=[];

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
 ok1=strcmp(ThisBac.birthtype, 'OK');    %birth ok
 ok2= strcmp(ThisBac.divtype, 'OK');    ;%division ok
 ok3=ThisBac.edges.edgesok;  %edges ok
 ok4=(length(frs)>0); 
 ok5=strcmp(S(ch).channels.RepClicks(j).fate, 'disassembled');
 ok6=S(ch).channels.AutoDivision(j).accepted;  %manual accept/reject
if ok1&ok2&ok3&ok4&ok5&ok6


%Select only mother cells, mother&first daughter, or all------------------
cond=0;
switch variant
    case '_all', cond=1;
    case '_mother',cond=(mod(log2(bacname),1)==0);
    case '_M1D1',  cond=((mod(log2(bacname),1)==0)|(mod(log2(bacname-1),1)==0));
end   
    
    firstframe=S(ch).channels.ReplicationCluster(j).PosKyTracCom.frames(1);
totalbaccount  =totalbaccount+1;
 %Select the fluorescent values-------------------------------------
%0) Select the time axis (normalized or not)
switch timeax
    case 'absolute', %time in minutes, starts at 0 per bacterium
        frames=minsperframe*S(ch).channels.ReplicationCluster(j).FluoPropsGen.frs_ext;
        frames=frames-frames(1);  
    case 'relative', %time in scaled units, 0 at init, 1 at ter
        frames=S(ch).channels.ReplicationCluster(j).Cycle.normtime_ext';
end
%------------------------------
%1) Select total count
all_labels_raw=S(ch).channels.ReplicationCluster(j).FluoPropsGen.signalcontent_ext;
%all_labels_raw=S(ch).channels.ReplicationCluster(j).Pos2DPreTrac.contentallspots';
%2) select  spots count (options available)
switch sptyp
case 1, peak_labels=S(ch).channels.ReplicationCluster(j).Pos2DPreTrac.contentallspots';
case 2, peak_labels=S(ch).channels.ReplicationCluster(j).FluoPropsGen.signalpeakval_ext*(2*pi*est.psf^2);
case 3, peak_labels=S(ch).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontent_ext;
case 4, peak_labels=S(ch).channels.ReplicationCluster(j).FluoPropsGen.fluospotscontentpadd_ext;
case 5, peak_labels=S(ch).channels.ReplicationCluster(j).FluoPropsGen.sumspotall_ext';
case 6, peak_labels=S(ch).channels.SecondFluoCluster(j).Pos2DPreTrac.contentallspots';
end


%In this section, data is corrected for the expected bleaching ------
abs_frames=S(ch).channels.ReplicationCluster(j).FluoPropsGen.frs_ext';  

%--For Tus-minus
peak_labels=Convert_fluorescence_to_labelcounts_DeltaTus(countsperlabel,peak_labels,abs_frames);
all_labels=Convert_fluorescence_to_labelcounts_DeltaTus(countsperlabel,all_labels_raw,abs_frames);

%peak_labels=Convert_fluorescence_to_labelcounts(countsperlabel,peak_labels,abs_frames);
%all_labels=Convert_fluorescence_to_labelcounts(countsperlabel,all_labels_raw,abs_frames);

%--------Converting to dimers
peak_labels=peak_labels/2;
all_labels=all_labels/2;


if 0
    plot(frames,peak_labels,'o-');
    [~]=ginput(1);   
end


%peak_labels=peak_labels/countsperlabel;
%all_labels= all_labels_raw/countsperlabel;
%------------------------------------------------------


%Jacob Kers, 22-1--2013----------------------------------------------
%In this section, we obtain an expected bacterium length. We use this to
%normalize the total count to a bacterium of average length. 
%This implies that if the concentration in the cytosol stays constant 
%and the spots content too, we expect a constant percentage
midlifebaclength=17;  %pixels; to calculate comparable expression levels
Mpars=ThisBac.edges.Mpars; %2nd order mid-position
Lpars=ThisBac.edges.Rpars;  %1st order right
Rpars=ThisBac.edges.Lpars;  %1st order left
frs=abs_frames;
left=(Mpars(1)*frs.^2+(Mpars(2)+Lpars(1))*frs+(Mpars(3)+Lpars(2)))';
right=(Mpars(1)*frs.^2+(Mpars(2)+Rpars(1))*frs+(Mpars(3)+Rpars(2)))';
baclength=(right-left)';
corfactor=baclength/midlifebaclength;
corfactor=1;
perc_labels=peak_labels./all_labels*100.*corfactor';  
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
%In this section, we concatenate bacterium data for later processing
lp=length(frames);
idxes=[1:lp]';
AllBax_times=[AllBax_times ; [idxes frames]];
AllBax_totals=[AllBax_totals ;[idxes all_labels']];
AllBax_spots=[AllBax_spots ; [idxes peak_labels']];
AllBax_perc=[AllBax_perc ; [idxes perc_labels']];
%AllBax_perc=[AllBax_perc ; [idxes baclength']];

%Collect replication and division times
tdiv=S(ch).channels.AutoDivision(j).divisiontime-...
     S(ch).channels.AutoDivision(j).birthtime;
trep=S(ch).channels.RepClicks(j).PosClick.lastframe-...
     S(ch).channels.RepClicks(j).PosClick.firstframe;
alldivtime(totalbaccount)=tdiv*minsperframe;
allreptime(totalbaccount)=trep*minsperframe; 
end
end
end

%-----------------------------------------------------------------------
%Run Binfillers: In this section, data is redivided over bins. 
%We choose to have one data point per bacterium per bin; Therefore the data is resampled.
%Contents of 'results'
% binresults.losig;          %one sigma below the average, per bin
% binresults.av;             % the average, per bin
% binresults.hisig;          %one sigma above the average, per bin
% binresults.binaxis;   %values of the non-empty bins, per bin
% binresults.binaxis_all; %value of  bin, per point
% binresults.binaxis_all_scat; %random number to spread bin points (for plotting purposes)
% binresults.binvalues_all; %values themselves
% binresults.bincounts=traceno;
binresultscell=Binfiller_General_IntPol(AllBax_times,AllBax_totals,binax,1);
binresultsspots=Binfiller_General_IntPol(AllBax_times,AllBax_spots,binax,1);
binresultsperc=Binfiller_General_IntPol(AllBax_times,AllBax_perc,binax,1);
%-----------------------------------------------------------------------


%---------------------------------------------------
%In this section, we make  stacked histograms of the data points over a period of
%normalized time, for example from t=0.2 to 0.8, where t=0 indicates
%initiation and t=1 indicates termination. The aim is to see if outliers
%affect our average & sigma too much. Note : this action only makes sense
%if we use relative times hence the condition
if strcmp(timeax,'relative');      %for cycle
slots=[-0.25:0.25:1.25];  %these are begin and end points of the relative time intervals
dum=Proc_Make_histograms_binresults(binresultsspots,slots);
end
%-------------------------------------------------------------


%-------------------------------------------------------------------------
%Here, we make Histograms of replication and division times and plot the
%result.
divhist=hist(alldivtime,absbinax);
rephist=hist(allreptime,absbinax);
figure; subplot(1,3,1);
bar(absbinax--minsperframe/2,divhist,0.5,'k'); hold on;
bar(absbinax+-minsperframe/2,rephist,0.5,'r'); hold on;
title(strcat(titlelabel,'-Replication and Division, N=',num2str(totalbaccount)));
xlabel('rel.time, minutes'); ylabel('counts');
axis([0 max(absbinax) 0 1.2*max([divhist rephist])]);
legend('division time', 'replication time');
%-------------------------------------------------


%General plot and save menu--------------------------------------------------------
%1a) Individual points and averages per bin--------------------------------------------------------------
%for plotting use, extract data-----------------
ax2D=binresultscell.binaxis_all;
axscat2D=binresultscell.binaxis_all_scat;
all_cell=binresultscell.binvalues_all;
all_spots=binresultsspots.binvalues_all;
all_perc=binresultsperc.binvalues_all;

loall=binresultscell.losig;
avall=binresultscell.av;
hiall=binresultscell.hisig;

lopeaks=binresultsspots.losig;
avpeaks=binresultsspots.av;
hipeaks=binresultsspots.hisig;

loperc=binresultsperc.losig;
avperc=binresultsperc.av;
hiperc=binresultsperc.hisig;
%----------------------------

subplot(1,3,2);
%1a)  total&peaks per bacterium per time point-----------------------------
plot(axscat2D,all_cell, 'bo', 'MarkerSize', 1.5); hold on;
plot(axscat2D,all_spots, 'ro', 'MarkerSize', 1.5); hold on;
title(strcat('Binned Fluorescence intensities-'));
xlabel('rel.time, minutes');
ylabel('label counts');
axis([minbintime max(binax) 0 7E5/countsperlabel]);
%1b totals peaks av-lo-hi
plot(binax,avall,'--ks','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',5); hold on;
plot(binax,loall,'k-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,hiall,'k-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
%1c)peaks av-lo-hi            
plot(binax,avpeaks,'--ks','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',5); hold on;
plot(binax,lopeaks,'k-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,hipeaks,'k-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
%--------------------------------------------------------------------------

%2 percentages------------------------------------------------------------
%2a percentage per bacterium per time point
subplot(1,3,3); 
plot(axscat2D,all_perc, 'ro', 'MarkerSize', 1.5); hold on;
title(strcat('Relative intensities-'));
xlabel('rel.time, minutes');
ylabel('spot intensities, % of total');

%2b perc av-lo-hi
plot(binax,avperc,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',5); hold on;
plot(binax,loperc,'r-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',5); hold on;
plot(binax,hiperc,'r-','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',5); hold on;
axis([minbintime max(binax) 0 100]);
%--------------------------------------------------------------------------




% %%------------------------HERE WE PLOT THE NUMBER OF MOLECULES (REL/ABS) IN CYTOPLASM
% %%------------------------AS WELL AS INSIDE THE FOCI

h=figure;
h1 = plot(axscat2D,all_cell, 'bo', 'MarkerSize', 1.5); hold on;
h1=h1(1);
h2 = plot(axscat2D,all_spots, 'ro', 'MarkerSize', 1.5); hold on;
h2=h2(1);
title(strcat('Increase of foci intensity after initiation , N=',num2str(totalbaccount)), 'fontsize', 16, 'fontweight', 'bold')
xlabel('Normalized Time (-)', 'fontsize', 16, 'fontweight', 'bold');
ylabel('Integrated intensity (-)', 'fontsize', 16, 'fontweight', 'bold');
%ylabel('Number of \beta_2 clamps (-)', 'fontsize', 16, 'fontweight', 'bold');
set(gca, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%axis([minbintime max(binax) 0 maxax]);
%1b totals peaks av-lo-hi
% 
% lmwrite('../week1/LeadDataTest.txt', [AlX AlY], 'delimiter', '\t', ...
% %          'precision', 3)


%Subtract ofset
% avpeaks=avpeaks-10;
% lopeaks=lopeaks-10;
% hipeaks=hipeaks-10;

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
%1c)peaks av-lo-hi            
h5 = plot(binax,avpeaks,'--ro','LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on;
h6 = plot(binax,lopeaks,'r-','LineWidth',3,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
plot(binax,hipeaks,'r-','LineWidth',3,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',5); hold on;
%axis([-0.2 0.2 0 4E5/countsperlabel]);
%avpeaks
%  dlmwrite('Foci_WT_DnaN_Avg.txt', [binax' avpeaks], 'delimiter', '\t','precision', 3)
%  dlmwrite('Foci_WT_DnaN_LOW.txt', [binax' lopeaks], 'delimiter', '\t','precision', 3)  
%  dlmwrite('Foci_WT_DnaN_HIGH.txt', [binax' hipeaks], 'delimiter', '\t','precision', 3) 
 
%   dlmwrite('Cyto_WT_DnaN_Avg.txt', [binax' avall], 'delimiter', '\t','precision', 3)
%  dlmwrite('Cyto_WT_DnaN_LOW.txt', [binax' loall], 'delimiter', '\t','precision', 3)  
%  dlmwrite('Cyto_WT_DnaN_HIGH.txt', [binax' hiall], 'delimiter', '\t','precision', 3) 
 

% dlmwrite('Foci_GammaMinus_DnaN_Avg.txt', [binax' avpeaks], 'delimiter', '\t','precision', 3)
% dlmwrite('Foci_GammaMinus_DnaN_LOW.txt', [binax' lopeaks], 'delimiter', '\t','precision', 3)  
% dlmwrite('Foci_GammaMinus_DnaN_HIGH.txt', [binax' hipeaks], 'delimiter', '\t','precision', 3) 

hleg1 = legend([h1 h2 h3 h4 h5 h6],{'Binned data points in cell', 'Binned data points in foci','Average intensity in cell', 'Standard deviation in cell', 'Average intensity in foci', 'Standard deviation in foci'});
set(hleg1,'Location','NorthWest')
      
lbl=strcat(initval.basepath,FileName_FociCyto);
%print(h, '-dpdf', '-r600',lbl)
%print(h, '-dpdf', '-r600','FociCyto_DnaN')




%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%-------Here we plot and write to file a selection of the absolute foci amount of beta_2 clamps data as
%%well as make a histogram and do a gaussian fit to it
N=length(axscat2D);
%xValues = zeros(N,1);
%yValues = zeros(N,1);
c=0;

for R = 1:N
    if (axscat2D(R) > 0.1) && (axscat2D(R) < 1) && (isnan(axscat2D(R)) == 0) && (isnan(all_spots(R)) == 0)
        c=c+1;
        xValuesFoci(c) = axscat2D(R);
        yValuesFoci(c) = all_spots(R);
    end
end
figure;
%subtract lower offset
yValuesFoci = yValuesFoci;
mean(yValuesFoci)
plot(xValuesFoci,yValuesFoci,'bo')


CombinedSeries = yValuesFoci;
N_CombinedSeries = length(CombinedSeries)
mean_Value = mean(CombinedSeries)
std_Value = std(CombinedSeries)

H =figure;
binWidth = 7; %This is the bin width
binCtrs = 0:binWidth:110; %Bin centers, depends on your data
%n=length(Series8(1,:));




counts = hist(CombinedSeries,binCtrs);

   %prob = counts / (n * binWidth);
   %prob = counts;
   prob = counts/trapz(binCtrs,counts); % here we normalize the histogram by dividing with the area under the curve. This will provide us with the PDF
   set(gca,'Color',[1,0.4,0.6])

   bar(binCtrs,prob);
   %set(gca, 'fontsize', 26, 'linewidth', 2, 'fontweight', 'bold');
   title('Signal in foci', 'fontsize', 18, 'fontweight', 'bold')
   ylabel('PDF (-)', 'fontsize', 26, 'fontweight', 'bold')
   xlabel('Number of \beta_2-mYPet molecules (-)', 'fontsize', 26, 'fontweight', 'bold')
   %a = annotation('textbox','Position',[0.65 0.65 0.35 0.16], 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('N = %d',N_CombinedSeries) );
   axis([0 100 0 4e-2])
   set(gca, 'XTick', [0 20 40 60 80 100]);
   %set(gca, 'YTick', [0 0.01 0.02 0.03 0.04]);
   set(gca,'TickLength',[0.02 0.02]);
   box off;
   set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
   
   %dlmwrite('Hist_Data_GammaMinus.txt', [binCtrs' prob'], 'delimiter', '\t','precision', 3)
   %dlmwrite('Hist_Data_WT.txt', [binCtrs' prob'], 'delimiter', '\t','precision', 3)
% 
% 
% meanValue = mean_Value;
% sigmaValue = std_Value;
%    std_error_Mean = sigmaValue./sqrt(N_CombinedSeries)
%    x=1:1:100;
%    g=(1/(sigmaValue.*sqrt(2*pi))).*2.7182818284.^(-0.5*(((x-meanValue)./sigmaValue).^2));%# pdf of the normal distribution
%    
%    PD = fitdist(CombinedSeries', 'normal')
%     hold on
%    plot(x,g,'k', 'linewidth', 3);
%    
%    
% a = annotation('textbox','Position',[0.65 0.65 0.45 0.05], 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('N = %d\nAvg Nr = %d\nStdDev = %d \%',N_CombinedSeries,round(meanValue),round(sigmaValue)) );
%    
% hleg = legend('Binned data','Gaussian fit');
%    set(hleg, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%    hold off 

  print(H, '-dpdf', '-r600','Histogram_dnaN_mYPet_AbsNr_Betas')

%    h_chi = chi2gof(counts)  
%    h_ktest = kstest(counts)
%    h_jbtest = jbtest(counts)%
%%---=---------------------------------------------------------------- 









%%-------------------HERE WE PLOT THE RATIO OF FLUORESCENCE IN CYTOPLASM
%%-------------------AND IN FOCI
h=figure;
% xTemp = binax2+scatax;
% yTemp = bincollector_perc_labels;
% h1 = plot(xTemp,yTemp, 'o', 'color' ,[0.4 0.4 0.4], 'MarkerSize', 2); 
h1 = plot(axscat2D,all_perc, 'o', 'color' ,[0.4 0.4 0.4], 'MarkerSize', 2); 
h1=h1(1);
hold on;
%title(strcat(titlelabel,'-Replication and Division, N=',num2str(totalbaccount)));
title(strcat('Ratio of fluorescence detected in foci\newlinewith respect to the whole cell, N=',num2str(totalbaccount)), 'fontsize', 16, 'fontweight', 'bold')
%xlabel('rel.time, minutes');
%ylabel('spot intensities, % of total');
xlabel('Normalized time (-)', 'fontsize', 16, 'fontweight', 'bold')
ylabel(['Fluorescence intensity ratio (%) \newline Foci / whole cell'], 'fontsize', 16, 'fontweight', 'bold')
set(gca, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');


%2b perc av-lo-hi
h2 = plot(binax,avperc,'--ro','LineWidth',3,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on;
            
h3 = plot(binax,loperc,'b-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerSize',4); hold on;
            
    plot(binax,hiperc,'b-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerSize',4); hold on;
axis([minbintime max(binax) 0 125]);

legend([h1 h2 h3],{'Binned data points', 'Average', 'Standard deviation'});

%a = annotation('textbox','Position',[0.65 0.65 0.45 0.16], 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('N = %d',N_CombinedSeries) );
%legend([h1 h2 h3 h4],{'p1','p2','p3','p4'}) ;

% hleg = legend('Data Points','Standard deviation','Average');
% % Make the text of the legend italic and color it brown
% set(hleg);


lbl=strcat(initval.basepath,FileName_Percentage);
%print(h, '-dpdf', '-r600',lbl)
%print(h, '-dpdf', '-r600','Percentage_DnaN')


% %%-----PLOTTING REPLICATION AND DIVISION TIMES SEPERATELY HISTOGRAM
% H1 =figure;
% binWidth =13; %This is the bin width
% binCtrs = 20:binWidth:165; %Bin centers, depends on your data
% 
% counts = hist(alldivtime,binCtrs);
% %h = histfit(CombinedSeries,7)
% 
%    %prob = counts / (n * binWidth);
%    %prob = counts;
%    %prob = counts/trapz(binCtrs,counts); % here we normalize the histogram by dividing with the area under the curve. This will provide us with the PDF
%    
%    bar(binCtrs,counts,'k');
% 
%    set(gca, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%    xlabel('Time (min)', 'fontsize', 16, 'fontweight', 'bold');
%    ylabel('Counts (-)', 'fontsize', 16, 'fontweight', 'bold');
%    title(strcat('Division time,   T =  {37} ^oC, N =  ',num2str(totalbaccount)), 'fontsize', 16, 'fontweight', 'bold')
%    AvgDivTime = round(mean(alldivtime));
%    a = annotation('textbox','Position',[0.65 0.65 0.85 0.06], 'fontsize', 14, 'linewidth', 1.5, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('t_{div} = %d min',AvgDivTime) );
% 
% print(H1, '-dpdf', '-r600','DivTimes_Hist_37degrees')
%    
% H2 =figure;
%    %binWidth =12;
%    counts = hist(allreptime,binCtrs);
%    %prob = counts/trapz(binCtrs,counts); % here we normalize the histogram by dividing with the area under the curve. This will provide us with the PDF
%    bar(binCtrs,counts,'r');
%    set(gca, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%    xlabel('Time (min)', 'fontsize', 16, 'fontweight', 'bold');
%    ylabel('Counts (-)', 'fontsize', 16, 'fontweight', 'bold');
%    title(strcat('Division and replication times,   T =  {37} ^oC, N =  ',num2str(totalbaccount)), 'fontsize', 16, 'fontweight', 'bold')
% 
% 
% AvgRepTime = round(mean(allreptime));
%  a = annotation('textbox','Position',[0.65 0.65 0.85 0.06], 'fontsize', 14, 'linewidth', 1.5, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('t_{rep} = %d min',AvgRepTime) );
% 
% print(H2, '-dpdf', '-r600','RepTimes_Hist_37Degrees')
%%-----PLOTTING REPLICATION AND DIVISION TIMES ON A SINGLE HISTOGRAM
H =figure;
% binWidth =20; %This is the bin width
% binCtrs = 65:binWidth:200; %Bin centers, depends on your data
 binWidth =10; %This is the bin width
binCtrs = 20:binWidth:150; %Bin centers, depends on your data

counts = hist(alldivtime,binCtrs);
%h = histfit(CombinedSeries,7)

   %prob = counts / (n * binWidth);
   %prob = counts;
   %prob = counts/trapz(binCtrs,counts); % here we normalize the histogram by dividing with the area under the curve. This will provide us with the PDF
   
   bar(binCtrs--minsperframe/2,counts,0.5,'k');
   hold on;
counts = hist(allreptime,binCtrs);
   %prob = counts/trapz(binCtrs,counts); % here we normalize the histogram by dividing with the area under the curve. This will provide us with the PDF
   bar(binCtrs-+minsperframe/2,counts,0.5,'r');
   hleg = legend('Division time', 'Replication time');
   set(hleg, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
   set(gca, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
   %xlabel('Time (min)', 'fontsize', 16, 'fontweight', 'bold');
   xlabel('Time (min)', 'fontsize', 26, 'fontweight', 'bold')
   %ylabel('Counts (-)', 'fontsize', 16, 'fontweight', 'bold');
   ylabel('Counts (-)', 'fontsize', 26, 'fontweight', 'bold')
   set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
   box off;
   set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
    set(gca,'TickLength',[0.02 0.02]);
   title(strcat('Division and replication times,   T =  {37} ^oC, N =  ',num2str(totalbaccount)), 'fontsize', 16, 'fontweight', 'bold')

AvgDivTime = round(mean(alldivtime));
AvgRepTime = round(mean(allreptime));
 a = annotation('textbox','Position',[0.65 0.65 0.85 0.06], 'fontsize', 14, 'linewidth', 1.5, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('t_{div} = %d min\nt_{rep} = %d min',AvgDivTime,AvgRepTime) );

lbl=strcat(initval.basepath,FileName_RepDivHistograms);

%print(H, '-dpdf', '-r600',lbl)
  
  
%print(H, '-dpdf', '-r600','RepDivHistograms_WT')




%Saving menu---------------------------------------------------------------
if actions.savedata
    %1) histograms
    histdata=[binax'-minsperframe/2 divhist' rephist'];
    filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',exp,variant,remark,'_DivRepHist.txt');
     dlmwrite(filname,histdata, 'delimiter', '\t');
     
     %2) binned averages
     avdata=[binax'-minsperframe/2 loall avall hiall lopeaks avpeaks hipeaks loperc avperc hiperc];
    filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',exp,variant,remark,'_SpotsAndTotals.txt');
     dlmwrite(filname,avdata, 'delimiter', '\t');
     
     %2) individual points per bin
     binpoints=[ax2D-minsperframe/2 axscat2D-minsperframe/2 all_cell all_spots all_perc];
    filname=strcat(initval.basepath,'Results_',initval.outname,'_Exp',exp,variant,remark,'_BinPoints.txt');
     dlmwrite(filname,binpoints, 'delimiter', '\t');
end
