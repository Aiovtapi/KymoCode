function [f1D,f2D]=Processing_stepwisetriplePeakFit(FL_spot,fitpresets,fluo_est,initval,plotit)
%Function first makes a double 1D fit, then 2D two Gaussians plus background
%psf is user-set and fixed

%output
% params1Dxy = [X0,X1,Y0,Y1, Background amplitude,Peak0, Peak1] from 2x1D XY estimate
% params2D = [X0,X1,Y0,Y1, Background amplitude,Peak0, Peak1] from 3x2D  estimate
%JacobKers, 2012
%------------------------------------------------
if nargin<5 %DEMO MODE setup
      close all
      pt='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\2013_08_14 FociEval\ImageDatabase\';
      switch 3  %choose demo picture        
          case 1, pth=strcat(pt,'SingleFocus\');
          case 2, pth=strcat(pt,'\NoFoci\');
          case 3, pth=strcat(pt,'\MultipleFoci\');
      end
    imno=10;
    imagenames=dir(strcat(pth ,'*.tif'));
    nm=strcat(pth,imagenames(imno).name);
    im0=double(imread(nm));  %load an image; including borders         
    [r,c]=size(im0);        %Crop to nonzero area
    [X,Y]=meshgrid(1:c,1:r);
    sel=find(im0~=0);
    lor=min(Y(sel));             hir=max(Y(sel)); 
    loc=min(X(sel));             hic=max(X(sel));
    im=im0(lor:hir,loc:hic); 
    
    %fit settings user%--------------
    initval.spotRoiHW=9;
    initval.skip2Danalysis=1;
    fitpresets.hwix=initval.spotRoiHW;   %spot-roi limits, vertical
    fitpresets.hwiy=9;  %spot-roi limits, lateral
    fitpresets.psf=1.3;  %user-set pointspread function
    
    %-------------------------
    plotit=0; 
    fluo_est=Processing_Fluorescence_PatternAnalysis(im); %Get some general pattern properties from this image
   [FL_spot,FL_cyto,FL_residu]=Processing_fluorescence_splitpic(im,fluo_est);
end

%---------------------------------------------------

[r,c]=size(FL_spot);
spotcurve=sum(FL_spot);
inspots=sum(spotcurve);

%1) precook image using initial estimates-------------------------------------------
%aim is to fit properly scaled but to preserve the number of counts in the
%spots
normim=(FL_spot)/inspots;  %normalized image supposedly only containing fluorescence counts!
prf_spot_nrm = spotcurve/inspots;
ax=(1:1:length(prf_spot_nrm));
[x,y]=meshgrid(1:c,1:r);


%1D-fit procedure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2 make initial estimates for normalized fit--------------------------
est.psf=fitpresets.psf;
est.N=0.5;  %peak content per peak 
est.x0=c/2-c/4;
est.x1=c/2+c/4;
est.y0=r/2;
est.y1=r/2;


%2x1D GaussFit on summed, normalized pixels!------------------------------
    [params1Dn,y1Fitn,y2Fitn, yAllFitn]=MLE_Two1D_Gaussian_RealDataFixPSF(ax,prf_spot_nrm,est);
    [y0,y1,y0pix,y1pix]=Get_perpendicular_Positions(params1Dn,FL_spot);
    params1Dxyn = [params1Dn(1),params1Dn(2),y0,y1,params1Dn(3), params1Dn(4)];
    %-------------------------------------------------------------------------
    %translate parameters and curves back to non-normalized units;-----------------------
    %in addition, peak values are translated back to integrated peak values
    %----------1 for xy 1D double gauss fit on 1D summed pixels:
    f1D.X0=params1Dxyn(1);  %X1 position, pixels
    f1D.X1=params1Dxyn(2);  %X2 position, pixels
    f1D.Y0=params1Dxyn(3);  %Y1 position, pixels
    f1D.Y1=params1Dxyn(4);  %Y2 position, pixels
    f1D.contentspot1=params1Dxyn(5)*inspots; %INTEGRATED counts spot 1
    f1D.contentspot2=params1Dxyn(6)*inspots;  %INTEGRATED counts spot 2
    f1D.contentallspots=inspots;  %checksum
    f1D.contentfitspots=(params1Dxyn(5)+params1Dxyn(6))*inspots;  %coverage
    

%2D fit Procedure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~initval.skip2Danalysis
    %1) update relevant estimates from 1D fit
    est2.x0=params1Dn(1); 
    est2.x1=params1Dn(2);
    est2.y0=y0pix;
    est2.y1=y1pix;
    est2.N0=params1Dn(3);
    est2.N1=params1Dn(4);
    %3---Do FIT!------------------------------
    if 0
    pcolor(normim); shading flat; colormap hot; hold on;
    plot(est2.x0,est2.y0,'bo'); hold on;
    plot(est2.x1,est2.y1,'go'); hold off;
    pause(0.2);
    end
    
    datafun = @(params)(sum (sum ((DoubleTwoDGauss_Fix1PSF(x,y,est.psf,params))))-sum (sum (normim.*log (DoubleTwoDGauss_Fix1PSF(x,y,est.psf,params))))); 
    options = optimset ('MaxFunEvals', 1E3, 'MaxIter', 1E4, 'TolFun', 1e-4);
    params2Dn= fminsearch (datafun, params1Dxyn, options);
    params2Dn(5)=abs(params2Dn(5));  %to ensure positive amplitudes (also in expectation expression)
    params2Dn(6)=abs(params2Dn(6));  %to ensure positive amplitudes (also in expectation expression)

    f2D.X0=params2Dn(1);  %X1 position, pixels
    f2D.X1=params2Dn(2);  %X2 position, pixels
    f2D.Y0=params2Dn(3);  %Y1 position, pixels
    f2D.Y1=params2Dn(4);  %Y2 position, pixels
    f2D.contentspot1=params2Dn(5)*inspots;  %INTEGRATED free label counts
    f2D.contentspot2=params2Dn(6)*inspots; %INTEGRATED counts spot 1
    f2D.contentallspots=inspots;  %checksum spots
    f2D.contentfitspots=(params2Dn(5)+params2Dn(6))*inspots;  %coverage
        
else
    f2D=f1D;
end

Plot_Intermediates(FL_spot,f1D,f2D,ax,spotcurve,y1Fitn,y2Fitn,yAllFitn,inspots, plotit);  %just plotting


function [y0,y1,y0pix,y1pix]=Get_perpendicular_Positions(params1Dn,im)
%2) make coarse estimates for perpendicular position----------
[r,c]=size(im);
x0_c=ceil(params1Dn(1));
x1_c=ceil(params1Dn(2));

if ((x0_c>1) & (x0_c<c))
prfy=im(:,x0_c);
[val,y0pix]=max(prfy);
else
prfy=mean(im')';
[val,y0pix]=max(prfy);
end
if y0pix==1, 
    y0pix=r/2; 
    y0=y0pix;
else
    y0= Refine_Peaks(prfy,y0pix, 0);
end

if ((x1_c>1) & (x1_c<c))
prfy=im(:,x1_c);
[val,y1pix]=max(prfy);
else
prfy=mean(im')';
[val,y1pix]=max(prfy);
end
if y1pix==1, 
    y1pix=r/2;
    y1=y1pix;
else
    y1= Refine_Peaks(prfy,y1pix, 0);
end

%----------------------------------------------------


    function dum=Plot_Intermediates(FLspot_narrow,f1D,f2D, ax,prf_all,y1Fitn,y2Fitn,yAllFitn,prf_spot_sum,show);
    %'pictures' ; %'pictures';% 'geometry'; 'pictures'; %'geometry'; %'pictures';%'geometry'; %'pictures';% 'geometry'  %'pictures'
    if show==1       
       %For showing fits
        %Y1=prf_spot_nrm*prf_spot_sum+prf_spot_sum;
        X1=ax;
        Y1 = prf_all;
        %+fluo_est.darklevel;
        y1Fit=y1Fitn*prf_spot_sum+min(prf_all);
        y2Fit=y2Fitn*prf_spot_sum+min(prf_all);%
        allFit=yAllFitn*prf_spot_sum+min(prf_all);

        subplot(2,1,1); %figure;
        plot(X1,Y1,'o-r','LineWidth',1.5);
        axis([X1(1) X1(end) 0 1.3*max(allFit)]);
        hold on;
        plot(X1,y1Fit,'b-*','LineWidth',1.5);
        plot(X1,y2Fit,'g-*','LineWidth',1.5);
        plot(X1,(allFit),'k-*','LineWidth',1.5);
        h = legend('Input Data','Fit of Data 1','Fit of Data 2','Fit1+Fit2',5);
        set(h,'FontSize',7);
        xlabel ('Position ','fontsize',14);
        ylabel ('Intensity [a.u.]','fontsize',14);
        hold off;
    
        % %plot pictures
        % pcolor(thischan_FL); shading flat; colormap hot; title('FL');
        % axis equal, 
        subplot(2,1,2); %figure; 
        pcolor(FLspot_narrow); shading flat; colormap hot;  title('Replicationcluster'); whitebg('white');   hold on;    
        plot(f1D.X0+0.5,f1D.Y0+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','b'); 
        plot(f1D.X1+0.5,f1D.Y1+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','g');
        plot(f2D.X0+0.5,f2D.Y0+0.5,'o','MarkerFaceColor', 'b', 'MarkerSize', 12, 'MarkerEdgeColor','b');
        plot(f2D.X1+0.5,f2D.Y1+0.5,'o','MarkerFaceColor', 'g', 'MarkerSize', 12, 'MarkerEdgeColor','g'); 
        %if spotno==2 
%         if (~isnan(f2D.Spot1OK) & ~isnan(f2D.Spot2OK))
%             plot([f2D.X0 f2D.X1]+0.5, [f2D.Y0 f2D.Y1]+0.5, '-');
%         end
         %[~]=ginput(1);
        pause(0.3);
        hold off;
        
        if show & 0
            hold on 
            outname=strcat(filename, int2str(j));
            saveas(gcf,outname,'jpg')
            F(k) = getframe;
            F(k) = getframe(gcf);
            image(F(k).cdata);
            colormap(F(k).colormap);
        end       
    end
    
    function [betterpeaks, betterpeaksvals]= Refine_Peaks(data,firstpeaks, plotit)
%This function refines initial estimates of peak positions via parabolic
%fitting of sub-sections around each initial estimate.
%input: 'data': 1D array of points
%input: 'first peaks': array of indices of first estimates (should be integer)
%output:'betterpeaks': array of refined positions (sub-unit resolution)
%output:'betterpeaksvals': array of refined peak values (parabola maxima)
%Jacob Kers, 2013
%--------------------------------------------------------
%Run settings-----------------
subsectionsize=5; %uneven!  %3 is minimal but often OK
%------------------------------------------------

betterpeaks=0*firstpeaks;
betterpeaksvals=0*firstpeaks;
lf=length(firstpeaks);
ld=length(data);

for i=1:lf % for all initial peak positions
    %1) pick a subsection around one peak pos; ensure borders to stay within
    %data limits
    mid=firstpeaks(i);
    lo=mid-(subsectionsize-1)/2;  lo=max(lo,1);
    hi=mid+(subsectionsize-1)/2; hi=min(hi,ld);
    indexaxis=[lo:1:hi]';
    subsection=data(indexaxis);
    %2) perform a parabolic fit and optionally plot results----------------
    subs2=(subsection-nanmean(subsection))/range(subsection); %scale
    ax2=indexaxis-mean(indexaxis);
    prms=polyfit(ax2,subs2,2); %fit
    prms(1)=prms(1)*range(subsection);
    prms(2)=prms(2)*range(subsection);
    prms(3)=prms(3)*range(subsection)+nanmean(subsection); %rescale
    peakpos=-prms(2)/(2*prms(1));                       %parabol max position
    peakval=prms(1)*peakpos^2+prms(2)*peakpos+prms(3);  %parabola max value
    betterpeaks(i)=peakpos+mean(indexaxis);
    betterpeaksvals(i)=peakval;
    %3 optional plot menu-------------------------------------------------
    if plotit  
    parabolfit=prms(1)*indexaxis.^2+prms(2)*indexaxis+prms(3);
    plot(indexaxis, parabolfit, '-r');
    end
end
