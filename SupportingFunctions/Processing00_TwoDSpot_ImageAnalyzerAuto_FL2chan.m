function Processing00_TwoDSpot_ImageAnalyzerAuto_FL2chan(exp,whichchan)
%Two-dim analysis of EcoliData for bacterial analysis; Jacob Kerssemakers,
%TNW-BN-ND lab 2012; developed for Charl Moolman; It is assumed there are 2
%channels available

if nargin<2, 
    exp='TEST';
    whichchan=1;
    plotit=1;
end

actions.loaddatabase=1; %default=1 (new analysis)

initval=A001_Images_Set_Experiment(exp);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------

[good, bad]=Processing_Measure_GoodDatabase(S)


%work through all replication cycles------------------------
[~,chan_no]=size(S);
difs=[]; doit=0;
%for i=1:chan_no  %for each channel
for i=1:chan_no  %for each channel
Rep=S(i).channels.ReplicationCluster;
Div=S(i).channels.AutoDivision;


ManRep=M(i).channels.RepClicks;
switch whichchan
    case 1  %take the first channel('DnaN')
    chanstk_FL=S(i).channels.chanstk_FL;
    kymo_FL=S(i).channels.kymo_FL;   
    case 2  %take the secondary fluorescence ('DnaX')
    chanstk_FL=S(i).channels.chanstk_FL2;
    kymo_FL=S(i).channels.kymo_FL2;
end
[~,repno]=size(Div)
dum=1;
for j=1:repno; %for each bacterium
close all
j
ThisRep=Rep(j);         %Current replication
ThisManRep=ManRep(j);
ThisDiv=Div(j);
frs=       ThisDiv.edges.frs; 
%check number of conditions
 ok1=strcmp(ThisDiv.birthtype, 'OK');    %birth ok
 ok2= strcmp(ThisDiv.divtype, 'OK');    ;%division ok
 ok3=ThisDiv.edges.edgesok;  %edges ok
 ok4=(length(frs)>0); 
 %ok5=ThisBac.accepted;
 ok6=strcmp(ThisManRep.fate, 'disassembled');
 ok9=ThisDiv.accepted;

 

if ok2&ok3&ok4&ok6&ok9
%Do extensive 1D and 2D spot analysis. result: %pre-fit, final fit: 
%[X0,X1,Y0,Y1, Background amplitude,Peak0, Peak1, spotno]
'bacteria to go' , good.birthdivedgesuser-doit
doit=doit+1;
[fluopropcurves,areasums,prefits,finalfits]=Processing_ClusterLife(i,j,ThisRep,ThisDiv,initval,chanstk_FL,whichchan); 
%--------------------------------------------------------------------------

%Update the -proper- database with analysis result
switch whichchan
    case 1
            %S(i).channels.ReplicationCluster(j).FluoPropsFin.
            S(i).channels.ReplicationCluster(j).FluoPropsFin.area_bac=fluopropcurves.area_bac;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.area_spot=fluopropcurves.area_spot;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.content_cytoplasm=fluopropcurves.content_cytoplasm;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.content_signal=fluopropcurves.content_signal;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.content_spots=fluopropcurves.content_spots;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.content_total=fluopropcurves.content_total;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.level_dark=fluopropcurves.level_dark;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.level_fluotreshold=fluopropcurves.level_fluotreshold;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.level_medianofmax=fluopropcurves.level_medianofmax;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.level_medianofsum=fluopropcurves.level_medianofsum;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.level_peak=fluopropcurves.level_peak;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.noise_dark=fluopropcurves.noise_dark;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.peak_xpos=fluopropcurves.peak_xpos;
            S(i).channels.ReplicationCluster(j).FluoPropsFin.peak_ypos=fluopropcurves.peak_ypos;

            S(i).channels.ReplicationCluster(j).AreaSumming.sumI1=areasums.I1;
            S(i).channels.ReplicationCluster(j).AreaSumming.sumI2=areasums.I2;
            S(i).channels.ReplicationCluster(j).AreaSumming.sumIall=areasums.Iall;

            % %Add extended position data JK:????????????????
            % S(i).channels.ReplicationCluster(j).PosKyTracCom.frames_ext=ThisRep.PosKyTracCom.frames_ext;  
            % S(i).channels.ReplicationCluster(j).PosKyTracCom.trackpos_ext=ThisRep.PosKyTracCom.frames_ext;  

            %Fluorecencent spot: Add area summation data
            S(i).channels.ReplicationCluster(j).FluoPropsGen.areasumI1=areasums.I1;
            S(i).channels.ReplicationCluster(j).FluoPropsGen.areasumI2=areasums.I2;
            S(i).channels.ReplicationCluster(j).FluoPropsGen.areasumIall=areasums.Iall;

            %Fluorecencent spots: add 2x 1D Gauss fit (along X) plus Y estimates
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.X0=prefits.X0;
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.X1=prefits.X1;
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.Y0=prefits.Y0;
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.Y1=prefits.Y1;
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.contentallspots=prefits.contentallspots;
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.contentspot1=prefits.contentspot1; %pixel contents spot 1
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.contentspot2=prefits.contentspot2; %pixel contents  spot 2
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.spot1OK=prefits.spot1OK;
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.spot2OK=prefits.spot2OK;

            %Fluorecencent spots: add 2D Gauss fit 
            S(i).channels.ReplicationCluster(j).Pos2DFinTrac.X0=finalfits.X0;
            S(i).channels.ReplicationCluster(j).Pos2DFinTrac.X1=finalfits.X1;
            S(i).channels.ReplicationCluster(j).Pos2DFinTrac.Y0=finalfits.Y0;
            S(i).channels.ReplicationCluster(j).Pos2DFinTrac.Y1=finalfits.Y1;
            S(i).channels.ReplicationCluster(j).Pos2DFinTrac.contentallspots=finalfits.contentallspots;%Background
            S(i).channels.ReplicationCluster(j).Pos2DFinTrac.contentspot1=finalfits.contentspot1; %pixel contents  spot 1
            S(i).channels.ReplicationCluster(j).Pos2DFinTrac.contentspot2=finalfits.contentspot2; %pixel contents  spot 2   
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.spot1OK=finalfits.spot1OK;
            S(i).channels.ReplicationCluster(j).Pos2DPreTrac.spot2OK=finalfits.spot2OK;
            %--------------c------------------------------------------------------------
    case 2
        %S(i).channels.SecondFluoCluster(j).FluoPropsFin.
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.area_bac=fluopropcurves.area_bac;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.area_spot=fluopropcurves.area_spot;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.content_cytoplasm=fluopropcurves.content_cytoplasm;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.content_signal=fluopropcurves.content_signal;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.content_spots=fluopropcurves.content_spots;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.content_total=fluopropcurves.content_total;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.level_dark=fluopropcurves.level_dark;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.level_fluotreshold=fluopropcurves.level_fluotreshold;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.level_medianofmax=fluopropcurves.level_medianofmax;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.level_medianofsum=fluopropcurves.level_medianofsum;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.level_peak=fluopropcurves.level_peak;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.noise_dark=fluopropcurves.noise_dark;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.peak_xpos=fluopropcurves.peak_xpos;
            S(i).channels.SecondFluoCluster(j).FluoPropsFin.peak_ypos=fluopropcurves.peak_ypos;

            S(i).channels.SecondFluoCluster(j).AreaSumming.sumI1=areasums.I1;
            S(i).channels.SecondFluoCluster(j).AreaSumming.sumI2=areasums.I2;
            S(i).channels.SecondFluoCluster(j).AreaSumming.sumIall=areasums.Iall;
         

            %Fluorecencent spots: add 2x 1D Gauss fit (along X) plus Y estimates
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.X0=prefits.X0;
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.X1=prefits.X1;
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.Y0=prefits.Y0;
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.Y1=prefits.Y1;
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.contentallspots=prefits.contentallspots;
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.contentspot1=prefits.contentspot1; %pixel contents spot 1
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.contentspot2=prefits.contentspot2; %pixel contents  spot 2
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.spot1OK=prefits.spot1OK;
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.spot2OK=prefits.spot2OK;

            %Fluorecencent spots: add 2D Gauss fit 
            S(i).channels.SecondFluoCluster(j).Pos2DFinTrac.X0=finalfits.X0;
            S(i).channels.SecondFluoCluster(j).Pos2DFinTrac.X1=finalfits.X1;
            S(i).channels.SecondFluoCluster(j).Pos2DFinTrac.Y0=finalfits.Y0;
            S(i).channels.SecondFluoCluster(j).Pos2DFinTrac.Y1=finalfits.Y1;
            S(i).channels.SecondFluoCluster(j).Pos2DFinTrac.contentallspots=finalfits.contentallspots;%Background
            S(i).channels.SecondFluoCluster(j).Pos2DFinTrac.contentspot1=finalfits.contentspot1; %pixel contents  spot 1
            S(i).channels.SecondFluoCluster(j).Pos2DFinTrac.contentspot2=finalfits.contentspot2; %pixel contents  spot 2   
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.spot1OK=finalfits.spot1OK;
            S(i).channels.SecondFluoCluster(j).Pos2DPreTrac.spot2OK=finalfits.spot2OK;
            %--------------c------------------------------------------------------------
end






%plot menu for diagnosis---------------------------------------------------
if initval.plotintermediateresults
X=ThisRep.PosKyTracCom.frames_ext;
YTot=fluopropcurves.content_signal;
Y1=fluopropcurves.content_spots;  %median excess
sel=find(Y1~=0); Y1=Y1(sel); X1=X(sel);

Y2=nansum([prefits.contentspot1.*prefits.spot1OK ; prefits.contentspot2.*prefits.spot2OK]);            %spot count 1D fit
sel=find(Y2~=0); Y2=Y2(sel); X2=X(sel); 

Y3=nansum([finalfits.contentspot1.*finalfits.spot1OK ; finalfits.contentspot2.*finalfits.spot2OK]);      %spot count 2D fits
sel=find(Y3~=0); Y3=Y3(sel); X3=X(sel);

Y4=areasums.Iall;                             %spot count, areafit
sel=find(Y4~=0); Y4=Y4(sel); X4=X(sel);

if S(i).channels.RepClicks(j).name>0  %no starters
mark=[0*max(Y2) 1.5*max(YTot)];
fr0=S(i).channels.RepClicks(j).PosClick.firstframe;
fr1=S(i).channels.RepClicks(j).PosClick.lastframe;

subplot(1,2,2);
plot(X,YTot,'b-o','LineWidth',4,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on
plot(X1,Y1,'r-o','LineWidth',4,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on 
plot(X2,Y2,'g-o','LineWidth',4,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','w',...
                'MarkerSize',8); hold on
plot(X3,Y3,'k-o','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on
plot(X4,Y4,'k-*','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',6); hold on
            
plot([fr0 fr0],mark, 'c-','LineWidth',4);
plot([fr1 fr1],mark, 'c-','LineWidth',4);
axis([X(1) X(end) mark]);

title('total intensity  and spot(s) intensity, 1D fits')
xlabel('Time (frames)','fontsize', 20, 'fontweight', 'bold');
ylabel ('Intensity (a.u.)','fontsize', 20, 'fontweight', 'bold');
set(gca, 'fontsize', 20, 'linewidth', 4, 'fontweight', 'bold','box','off');
set(gca,'TickLength',[0.02 0.02]);

%legend('Total Intensity', 'median excess', 'area summing' ,'1D Gauss','2D Gauss');
%legend('Total Intensity', 'Median excess', 'Area summing' ,'1D Gauss');
hleg = legend('Total Intensity', 'Median excess','1D2XGauss','2D2XGauss', 'AreaSum' );
set(hleg,'FontSize',8, 'linewidth', 0, 'fontweight', 'bold','location','NorthWest');


 dum=Processing_Map_FluorescenceAutoDiv(ThisDiv,ThisRep,chanstk_FL,kymo_FL,initval);

end 
 
if ~initval.skippicturesavingbecauseCharlislogginghislive 
    FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'Double1D_GaussianFitting/FitingSummary/'));
    if FolderExistence == 0
        mkdir(strcat(initval.basepath,initval.FiguresFolder,'Double1D_GaussianFitting/FitingSummary/'));
    end
    FileNameToSaveFittingResults = strcat('ChNo',num2str(i),'_BacNo',num2str(j),'_FittingResultSummary');
    h=gcf;
    print(h, '-dpng', '-r600',strcat(initval.basepath,initval.FiguresFolder,'Double1D_GaussianFitting/FitingSummary/',FileNameToSaveFittingResults));
end
%------------------------------------------------------
end
else
    [~,repno]=size(Rep)
    switch whichchan        
        case 1, S(i).channels.ReplicationCluster(j).dum=[];
        case 2, S(i).channels.SecondFluoCluster(j).dum=[];
    end
end
end
end
 
save(outname, 'S', '-append');
disp('done');


