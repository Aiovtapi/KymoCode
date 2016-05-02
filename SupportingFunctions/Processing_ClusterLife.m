function [fluopropcurves,areasums,prefits,finalfits]=Processing_ClusterLife(ChNo,BacNo,ThisRep,ThisDiv,initval,chanstk_FL,chantype,chanstk_BF);
%pre-fit, final fit, [X0,X1,Y0,Y1, Background amplitude,Peak0, Peak1]

%-------demo settings-----
show=0;
plotit=show;
%mov=1*show;
mov=0;
savebacpics=1;
%--------------------------

%fit settings user%--------------
fitpresets.hwix=initval.spotRoiHW;   %spot-roi limits, vertical
fitpresets.hwiy=initval.spotregionhalfwidth;  %spot-roi limits, lateral
fitpresets.psf=initval.psf;  %user-set pointspread function 
%-------------------------

%lbl=num2str(ThisRep.name);

%Define a fixed-width area around the center-of-mass tracked position-----
%frs=ThisRep.FluoPropsGen.frs_ext;

%RDL -- for prolonged time bacpics -- frs2 is size of prolonged (Nframes = 10) edges fit
frs=ThisDiv.edges.frs;

% left=ThisRep.PosKyTracCom.trackpos_ext-fitpresets.hwi; 
% right=ThisRep.PosKyTracCom.trackpos_ext+fitpresets.hwi ; 


%Find edges based on parameters fitted results of the edges

Mpars=ThisDiv.edges.Mpars; %2nd order mid-position
Lpars=ThisDiv.edges.Rpars;  %1st order right
Rpars=ThisDiv.edges.Lpars;  %1st order left

%%I ASSUME THESE EDGES DEFINE THE EFFECT ROI WHERE THE PROGRAM LOOKS TO FIT
left=(Mpars(1)*frs.^2+(Mpars(2)+Lpars(1))*frs+(Mpars(3)+Lpars(2)))';
right=(Mpars(1)*frs.^2+(Mpars(2)+Rpars(1))*frs+(Mpars(3)+Rpars(2)))';
spotpos=ThisRep.PosKyTracCom.trackpos_ext;
spotpos=[spotpos ;zeros(length(frs)-length(spotpos),1)];



%define some more arrays for upcoming analysis
lfr=length(frs); 
% nofits=zeros(lfr,3); 
% prefits=zeros(lfr,7);  %parameters of 2x 1D fit per frame
% finalfits=zeros(lfr,7); %parameters of 2D fit per frame
F0_area_bac=zeros(lfr,1);
F0_area_spot=zeros(lfr,1);
F0_content_cytoplasm=zeros(lfr,1);
F0_content_signal=zeros(lfr,1);
F0_content_spots=zeros(lfr,1);
F0_content_total=zeros(lfr,1);
F0_level_dark=zeros(lfr,1);
F0_level_fluotreshold=zeros(lfr,1);
F0_level_medianofmax=zeros(lfr,1);
F0_level_medianofsum=zeros(lfr,1);
F0_level_peak=zeros(lfr,1);
F0_noise_dark=zeros(lfr,1);
F0_peak_xpos=zeros(lfr,1);
F0_peak_ypos=zeros(lfr,1);
F1X0=zeros(lfr,1);
F1X1=zeros(lfr,1);
F1Y0=zeros(lfr,1);
F1Y1=zeros(lfr,1);
F1A=zeros(lfr,1);
F1S1=zeros(lfr,1);
F1S2=zeros(lfr,1);
F1OK1=zeros(lfr,1);
F1OK2=zeros(lfr,1);
F2X0=zeros(lfr,1);
F2X1=zeros(lfr,1);
F2Y0=zeros(lfr,1);
F2Y1=zeros(lfr,1);
F2A=zeros(lfr,1);
F2S1=zeros(lfr,1);
F2S2=zeros(lfr,1);
F2OK1=zeros(lfr,1);
F2OK2=zeros(lfr,1);
%----------------------------------------


%%-----------------------------------------------------------------
% Run a loop on the bacterial frames & positions to fetch
%corresponding images %later to be expand to twoD fit)
dif1=[]; dif2=[]; 
% 

% Lengthen the bacterial view by certain amount of pixels:

AdditionalPixels=0;

  right=right+AdditionalPixels;
   left=left+AdditionalPixels;

for k=1:lfr  %for all frames
lfr-k;
if show==1,close all; end

%get input for 'this' frame --------------------------------
% fitpresets.level_dark=ThisRep.FluoPropsGen.level_dark_ext(k);
% fitpresets.freelabellevel=ThisRep.FluoPropsGen.fluospinelevel_ext(k);

%border check-------------------
thisfr=max([1 round(frs(k))]);
thisleft=round(max([(left(k)) 1]));    
thisright=round(min([(right(k)) initval.kymolength]));

if thisfr<=size(chanstk_FL,3)
thischan_FL=fliplr(squeeze(chanstk_FL(:,:,thisfr)));
else
thischan_FL=zeros(size(chanstk_FL,1),size(chanstk_FL,2));
end

if thisfr<=size(chanstk_BF,3)
thischan_BF=fliplr(squeeze(chanstk_BF(:,:,thisfr)));
else
thischan_BF=zeros(size(chanstk_BF,1),size(chanstk_BF,2));
end


% RDL -- Assumed that spotpos remains constant in
% extended area when Bacpics are prolonged

Du=nonzeros(spotpos);
if spotpos(k)==0
thisspot=Du(end);
else 
thisspot=spotpos(k);
end

lf_nr=round(max([thisspot-fitpresets.hwix 1])); 
rht_nr=round(min([thisspot+fitpresets.hwix initval.kymolength]));     

border_ok=(thisright>thisleft+6)&&(rht_nr>lf_nr+6);       


if 0
    %%plots a single frame of the current channel that is being analysed
    figure(1);
    imagesc(thischan_FL)
    axis equal tight
    box on
    colormap(hot)
    [~]=ginput(1);
    %--------------
end

%create ROI;%---------------------Do 1D and 2D analysis------------- --------------
if border_ok==1  %if this is a reasonable ROI to fit
            
    % RDL -- making the bac broader to observe dif site (mostly at pole)
    % correctly in oriZ-dif strain.

        baclengths=right-left;
        maxbaclength=max(baclengths);
    
    
        lo=int32(initval.kymohwidth-fitpresets.hwiy+1); 
        hi=int32(initval.kymohwidth+fitpresets.hwiy+1); 
        lf=int32(thisleft);
        rht=int32(thisright);

        
        FL_narrow=thischan_FL(initval.kymohwidth-fitpresets.hwiy+1:initval.kymohwidth+fitpresets.hwiy+1,lf_nr:rht_nr); %this is a narrow region around the spot
        [Ys,Xs]=size(thischan_FL);
        % make sure 
        if lo<0 || hi>Ys || lf<0 || rht>Xs
            FL_full=zeros(hi-lo,rht-lf);
            BF_full=zeros(hi-lo,rht-lf);
        else
        FL_full=thischan_FL(lo:hi,lf:rht); %This is the full bacterial area;  
        BF_full=thischan_BF(lo:hi,lf:rht);
        end
        
       if savebacpics
           SavePaddedImage(FL_full,maxbaclength,ChNo,BacNo,k,initval,chantype,thischan_FL,BF_full);
       end       
        %Get some general pattern properties from the image; we use 
        %a restricted area around the replication cluster
        fluo_est=Processing_Fluorescence_PatternAnalysis(FL_narrow); 
        fluofull=Processing_Fluorescence_PatternAnalysis(FL_full);  
        fluofull=orderfields(fluofull);
        %Prepare for spotfitting: based on the fluorescence estimate
        %levels, subtract contributions from the cytoplasm
        [FL_spot,~]=Processing_fluorescence_splitpic(FL_narrow,fluo_est);  
        simprops.psf=1.5;
        [r,c]=size(FL_spot);
        
        %Add an insignificant spot to stabilize fit....
        %FL_fakespot=B201_Bacsimfocus(FL_spot, 10000, 3, 3, simprops,0);  
        
        
        
        %Do 1D and 2D analysis; evaluate results      
        [f1D,f2D]=Processing_stepwisetriplePeakFit(FL_spot,fitpresets,fluo_est,initval,plotit); 
        %pre-fit, finalfit ; note: we do this on a narrow, fixed width
        %area; therfore, the 'fluoprops' are not stored        
        
        f1D=Processing_EvalGaussFits(f1D,FL_spot,fluofull,initval);  %Evaluate results
        f2D=Processing_EvalGaussFits(f2D,FL_spot,fluofull,initval);
        
        [I1, I2,I_all]=Get_spotarea_Intensity(f1D,FL_spot,initval);
        %Get Intensities per spot and total (overlap corrected); 
        
        %Calculate back to absolute 'kymograph 'x-coordinates
        if 0  %if full area was used for fit
        f1D.X0=f1D.X0-1+thisleft;  %transform to absolute coordinates
        f1D.X1=f1D.X1-1+thisleft;  %transform to absolute coordinates
        f2D.X0=f2D.X0-1+thisleft;  %transform to absolute coordinates
        f2D.X1=f2D.X1-1+thisleft;  %transform to absolute coordinates
        else %if narrow area was used for fit
        f1D.X0=f1D.X0-1+(thisspot-fitpresets.hwix);  %transform to absolute coordinates
        f1D.X1=f1D.X1-1+(thisspot-fitpresets.hwix);  %transform to absolute coordinates
        f2D.X0=f2D.X0-1+(thisspot-fitpresets.hwix);  %transform to absolute coordinates
        f2D.X1=f2D.X1-1+(thisspot-fitpresets.hwix);  %transform to absolute coordinates   
        end
        
else   
[fluofull,f1D,f2D,I1, I2,I_all]=Set_as_bad_point;
end

%collect data points---------
F0_area_bac(k)=fluofull.area_bac;
F0_area_spot(k)=fluofull.area_spot;
F0_content_cytoplasm(k)=fluofull.content_cytoplasm;
F0_content_signal(k)=fluofull.content_signal;
F0_content_spots(k)=fluofull.content_spots;
F0_content_total(k)=fluofull.content_total;
F0_level_dark(k)=fluofull.level_dark;
F0_level_fluotreshold(k)=fluofull.level_fluotreshold;
F0_level_medianofmax(k)=fluofull.level_medianofmax;
F0_level_medianofsum(k)=fluofull.level_medianofsum;
F0_level_peak(k)=fluofull.level_peak;
F0_noise_dark(k)=fluofull.noise_dark;
F0_peak_xpos(k)=fluofull.peak_xpos;
F0_peak_ypos(k)=fluofull.peak_ypos;

F1X0(k)=f1D.X0;
F1X1(k)=f1D.X1;
F1Y0(k)=f1D.Y0;
F1Y1(k)=f1D.Y1;
F1A(k)=f1D.contentallspots;
F1S1(k)=f1D.contentspot1;
F1S2(k)=f1D.contentspot2;
F1OK1(k)=f1D.Spot1OK;
F1OK2(k)=f1D.Spot2OK;

F2X0(k)=f2D.X0;
F2X1(k)=f2D.X1;
F2Y0(k)=f2D.Y0;
F2Y1(k)=f2D.Y1;
F2A(k)=f2D.contentallspots;
F2S1(k)=f2D.contentspot1;
F2S2(k)=f2D.contentspot2;
F2OK1(k)=f2D.Spot1OK;
F2OK2(k)=f2D.Spot2OK;

nofits(k,:)=[I1,I2,I_all];  %estimated 'no-fit' based on area estimates via f1

%--------------
%'pictures' ; %'pictures';% 'geometry'; 'pictures'; %'geometry'; %'pictures';%'geometry'; %'pictures';% 'geometry'  %'pictures'
if ((show==1) & (exist('FLspot')==1)), wht='pictures'; else wht=0;end
switch wht
case 'pictures'
% %plot pictures
% pcolor(thischan_FL); shading flat; colormap hot; title('FL');
% axis equal,
if (exist('FLspot')==1)
whitebg('white');
subplot(2,1,1); 
lbl=num2str(k);
title(strcat('ReplicationCluster\_',lbl),'fontsize', 16 );
set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold','box','off');
ylim([0 15e4]);
subplot(2,1,2); %figure; 

pcolor(FL_spot); shading flat; colormap hot; title('FL');
set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
xlabel ('Position (px) ','fontsize', 26, 'fontweight', 'bold');
ylabel ('Position (px)','fontsize', 26, 'fontweight', 'bold');
%axis([1 2*fitpresets.hwi+1 1 2*fitpresets.hwiy+1]), 
title(strcat('FrameNr\_',num2str(k)),'fontsize', 16);hold on
whitebg('white');
%if spotno==2 
if (~isnan(f2D.Spot1OK) & ~isnan(f2D.Spot2OK))
%plot([f2D.X0 f2D.X1]+0.5, [f2D.Y0 f2D.Y1]+0.5, '-');
end

FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'Double1D_GaussianFitting/'));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'Double1D_GaussianFitting/'));
end
FileNameToSaveFittingResults = strcat('ChNo',num2str(ChNo),'_BacNo',num2str(BacNo),'_ReplicationCluster',lbl,'_FrameNr',num2str(k));

h=gcf;
print(h, '-dpng', '-r450',strcat(initval.basepath,initval.FiguresFolder,'Double1D_GaussianFitting/',FileNameToSaveFittingResults));
end
hold off;


if mov
F(k) = getframe(gcf);
image(F(k).cdata);
colormap(F(k).colormap);
end
dum=1;


case 'geometry'
plot(dif,'-o'); hold on; pause(0.1);
%axis([0 1 -90 90]);
title('pattern geometry');
% xlabel ('eccentricity');
% ylabel('angle, degrees');
ylabel ('eccentricity');
xlabel('time,frames');
[~]=ginput(1);
end
end


%Finally, fill output fields with parameters
fluopropcurves.area_bac=F0_area_bac;
fluopropcurves.area_spot=F0_area_spot;
fluopropcurves.content_cytoplasm=F0_content_cytoplasm;
fluopropcurves.content_signal=F0_content_signal;
fluopropcurves.content_spots=F0_content_spots;
fluopropcurves.content_total=F0_content_total;
fluopropcurves.level_dark=F0_level_dark;
fluopropcurves.level_fluotreshold=F0_level_fluotreshold;
fluopropcurves.level_medianofmax=F0_level_medianofmax;
fluopropcurves.level_medianofsum=F0_level_medianofsum;
fluopropcurves.level_peak=F0_level_peak;
fluopropcurves.noise_dark=F0_noise_dark;
fluopropcurves.peak_xpos=F0_peak_xpos;
fluopropcurves.peak_ypos=F0_peak_ypos;

prefits.X0=F1X0;
prefits.X1=F1X1;
prefits.Y0=F1Y0;
prefits.Y1=F1Y1;
prefits.contentspot1=F1S1;
prefits.contentspot2=F1S2;
prefits.contentallspots=F1A;
prefits.spot1OK=F1OK1;
prefits.spot2OK=F1OK2;

finalfits.X0=F2X0;
finalfits.X1=F2X1;
finalfits.Y0=F2Y0;
finalfits.Y1=F2Y1;
finalfits.contentspot1=F2S1;
finalfits.contentspot2=F2S2;
finalfits.contentallspots=F2A;
finalfits.spot1OK=F2OK1;
finalfits.spot2OK=F2OK2;

%finally, store in cluster
areasums=[];
areasums.I1=nofits(:,1);
areasums.I2=nofits(:,2);
areasums.Iall=nofits(:,3);


%Test plot menu
if initval.plotintermediateresults
    
    %%%Could THIS BE IT?>
    close(gcf);
    scrsz = get(0,'ScreenSize');
    figure('Position',[100 100 scrsz(3)/1.3 scrsz(4)/1.3]);
    subplot(1,2,1);
    plot(frs,left, 'c-','LineWidth',4); hold on;
    plot(frs,right, 'c-','LineWidth',4);
    h1 = plot(frs, spotpos, 'k-o','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','w',...
                'MarkerSize',6);
    h2 = plot(frs,prefits.X0.*prefits.spot1OK, 'b-o','LineWidth',4,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','w',...
                'MarkerSize',6);
    h3 = plot(frs,prefits.X1.*prefits.spot2OK, 'r-o','LineWidth',4,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',6);
    set(gca, 'fontsize', 20, 'linewidth', 4, 'fontweight', 'bold','box','off');
    set(gca,'TickLength',[0.02 0.02]);
    xlabel('Time (frames)','fontsize', 20, 'fontweight', 'bold');
    ylabel('Postion (px)','fontsize', 20, 'fontweight', 'bold');
    title('Positions different methods', 'fontsize', 14);
    
    hleg = legend([h1 h2 h3],{'SpotPos', 'PreFitsX0','PreFitsX1',});
    set(hleg,'FontSize',8, 'linewidth', 0, 'fontweight', 'bold','location','NorthWest'); 
    pause(0.01);
end
%---------------------


dum=1;
if mov
close all;
movie(F,1,3);
outname=strcat(outpath,outchan,'Replicationcluster',lbl,'_mov');
movie2avi(F,outname);
end
dum=1;
end


function [I1, I2,I_all, N1,N2,Nall]=Get_spotarea_Intensity(f1,FLspot,initval);
%This function obtains the summed intensity around predetermined maxima
%points; %prefits: x0,x1,y0,y1, b0,N0, N1
%first peak is brightest
%FLspot: image with spots between edges of bacterium


hw=initval.sumHW; %area width is (hw*2+1);

if 0 %plot menu--------------------------------------------------
close all
pcolor(FLspot); shading flat;

pcolor(FLspot); shading flat; colormap hot; title('FL');hold on
whitebg('white');
plot(f1.X0+1+0.5,f1.Y0+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','b');
plot(f1.X1+1+0.5,f1.Y1+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','g');
[~]=ginput(1);
end

[r,c]=size(FLspot);

%border patrol
x1=max(f1.X0+1, 1+hw); x1=min(x1, c-hw);
y1=max(f1.Y0, 1+hw); y1=min(y1, r-hw);
x2=max(f1.X1+1, 1+hw); x2=min(x2, c-hw);
y2=max(f1.Y1, 1+hw); y2=min(y2, r-hw);

%area coordinates (might overlap)
w1=round([x1-hw:1:x1+hw]);
h1=round([y1-hw:1:y1+hw]);
w2=round([x2-hw:1:x2+hw]);
h2=round([y2-hw:1:y2+hw]);

%find common points
w1w2=[w1 w2]; [~, mw, ~] = unique(w1w2); cmmw=w1w2(mw);
ovl_lo=length(cmmw)-(2*hw+1)+1; wcommons=w1(ovl_lo:2*hw+1);
h1h2=[h1 h2]; [~, mh, ~] = unique(h1h2); cmmh=h1h2(mh);
ovl_lo=length(cmmh)-(2*hw+1)+1; hcommons=h1(ovl_lo:2*hw+1);

%areas
    area1=FLspot(h1,w1);
    area2=FLspot(h2,w2);
    area3=FLspot(hcommons,wcommons);  %overlapping area!

%pixel numbers, this is later needed for background subtraction
N1=(2*hw+1)^2; N2=N1;
N3=length(area3(:));
Nall=N1+N2-N3;

%summed intensities; 
I1=sum(area1(:));
I2=sum(area2(:));
I3=sum(area3(:));

if ~isempty(I3)
I_all=I1+I2-I3;
else
I_all=I1+I2;
end

end

function [fluofull,f1D,f2D,I1, I2,I_all]=Set_as_bad_point;
%Dummy function to set NaN vals
fluofull.area_bac=NaN;
fluofull.area_spot=NaN;
fluofull.content_cytoplasm=NaN;
fluofull.content_signal=NaN;
fluofull.content_spots=NaN;
fluofull.content_total=NaN;
fluofull.level_dark=NaN;
fluofull.level_fluotreshold=NaN;
fluofull.level_medianofmax=NaN;
fluofull.level_medianofsum=NaN;
fluofull.level_peak=NaN;
fluofull.noise_dark=NaN;
fluofull.peak_xpos=NaN;
fluofull.peak_ypos=NaN;

f1D.X0=NaN;
f1D.X1=NaN;
f1D.Y0=NaN;
f1D.Y1=NaN;
f1D.contentallspots=NaN;
f1D.contentspot1=NaN;
f1D.contentspot2=NaN;
f1D.Spot1OK= NaN;
f1D.Spot2OK= NaN;

f2D.X0=NaN;
f2D.X1=NaN;
f2D.Y0=NaN;
f2D.Y1=NaN;
f2D.contentallspots=NaN;
f2D.contentspot1=NaN;
f2D.contentspot2=NaN;
f2D.Spot1OK= NaN;
f2D.Spot2OK= NaN;
f2D.Spot1OK= NaN;
f2D.Spot2OK= NaN;

I1=NaN;
I2=NaN;
I_all=NaN;

fluo.level_dark=NaN;
fluo.signalcontent=NaN;
fluo.maxspinelevel=NaN;
fluo.level_medianofsum=NaN;
fluo.spotscontent=NaN;
fluo.peakval=NaN;
fluo.peaky=NaN;
fluo.peakx=NaN;
end


function SavePaddedImage(FL_full,maxbaclength,ChNo,BacNo,k,initval,chantype,ChanFL,BF_full)
       
       % FLUORESCENCE IMAGE       
       %padd to equal lengths
       mx=uint32(ceil(maxbaclength))+1;
       [rs,cs]=size(FL_full);
       padd=ceil((mx-cs)/2);
       padstrip=zeros(rs,padd); 
       Adu=[padstrip FL_full padstrip];
       if size(Adu,2) == mx
       Savim=[padstrip FL_full padstrip];
       else 
       Adi=size(Adu,2)-mx;
       padstrip2=zeros(rs,padd-Adi);
       Savim=[padstrip FL_full padstrip2];
       end 
       [rs1,cs1]=size(Savim);
       padd2=ceil((cs1-rs1)/2);
       padstripupdown=zeros(padd2,cs1);
       Adum=[padstripupdown; Savim; padstripupdown];
       if size(Adum,1)~=size(Adum,2)
       Adif=size(Adum,1)-size(Adum,2);
       padstripupdown2=zeros(padd2-Adif,cs1);
       Savim=[padstripupdown; Savim; padstripupdown2];
       else
       Savim=[padstripupdown; Savim; padstripupdown];
       end
       
       % BRIGHTFIELD IMAGE
       %padd to equal lengths
%        mx=uint32(ceil(maxbaclength))+1;
%        [rs,cs]=size(BF_full);
%        padd=ceil((mx-cs)/2);
%        padstrip=zeros(rs,padd); 
%        Adu=[padstrip BF_full padstrip];
%        if size(Adu,2) == mx
%        Savim_BF=[padstrip BF_full padstrip];
%        else 
%        Adi=size(Adu,2)-mx;
%        padstrip2=zeros(rs,padd-Adi);
%        Savim_BF=[padstrip BF_full padstrip2];
%        end 
%        [rs1,cs1]=size(Savim_BF);
%        padd2=ceil((cs1-rs1)/2);
%        padstripupdown=zeros(padd2,cs1);
%        Adum=[padstripupdown; Savim_BF; padstripupdown];
%        if size(Adum,1)~=size(Adum,2)
%        Adif=size(Adum,1)-size(Adum,2);
%        padstripupdown2=zeros(padd2-Adif,cs1);
%        Savim_BF=[padstripupdown; Savim_BF; padstripupdown2];
%        else
%        Savim_BF=[padstripupdown; Savim_BF; padstripupdown];
%        end
       
       %Savim=Savim(:,1:mx);  %just to make sure
       Chl=num2str(ChNo,'%02.0f');
       Baxx=num2str(BacNo,'%04.0f');
       Imm=num2str(k,'%03.0f');
       baclabel= strcat('Fluo', num2str(chantype),....
            'Chan', num2str(Chl),...
            'Bac', num2str(Baxx));
%        baclabel_BF= strcat('BriFi', num2str(chantype),....
%             'Chan', num2str(Chl),...
%             'Bac', num2str(Baxx));
       bacpth=strcat(initval.basepath,initval.FiguresFolder,...
            'BacPics/',initval.viewchannel,'/',...
            baclabel,...
            '/');
        
%        bacpth_BF=strcat(initval.basepath,initval.FiguresFolder,...
%             'BacPics/Brightfield/',...
%             baclabel,...
%             '/'); 
        
       if ~isdir(bacpth), mkdir(bacpth); end
%        if ~isdir(bacpth_BF), mkdir(bacpth_BF); end
                   %'PadX',num2str(PadsizeX),...
            %'PadY',num2str(PadsizeY),...
            
        lbl_Fig=strcat(bacpth,...
            baclabel,...         
            'Im',num2str(Imm),...
            '.tif'); %bacterium picture
               
%         lbl_Fig_BF=strcat(bacpth_BF,...
%             baclabel_BF,...         
%             'Im',num2str(Imm),...
%             '.tif'); %bacterium picture
        
        bacim= uint16(Savim); 
        %bacimBF= uint16(Savim_BF);
        imwrite(bacim,lbl_Fig,'tif');
        %imwrite(bacimBF,lbl_Fig_BF,'tif');
end
