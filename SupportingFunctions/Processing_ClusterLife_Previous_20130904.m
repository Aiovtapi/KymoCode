function [fluoprops,areasums,prefits,finalfits]=Processing_ClusterLife(ThisRep,ThisDiv,initval,chanstk_FL);
%pre-fit, final fit, [X0,X1,Y0,Y1, Background amplitude,Peak0, Peak1]

%-------demo settings-----
show=1;
plotit=show;
%mov=1*show;
mov=0;
%--------------------------

areasums=1;

%administration-----------------
outpath=initval.basepath;
outchan=initval.outname; 

%fit settings user%--------------
fitpresets.hwi=initval.spotRoiHW;   %spot-roi limits, vertical
fitpresets.hwiy=9;  %spot-roi limits, lateral
fitpresets.psf=1.4;  %user-set pointspread function 
%-------------------------

%lbl=num2str(ThisRep.name);

%Define a fixed-width area around the center-of-mass tracked position-----
frs=ThisRep.PosKyTracCom.frames_ext; 
% left=ThisRep.PosKyTracCom.trackpos_ext-fitpresets.hwi; 
% right=ThisRep.PosKyTracCom.trackpos_ext+fitpresets.hwi ; 


%Find edges based on parameters fitted results of the edges
Mpars=ThisDiv.edges.Mpars; %2nd order mid-position
Lpars=ThisDiv.edges.Rpars;  %1st order right
Rpars=ThisDiv.edges.Lpars;  %1st order left

left=(Mpars(1)*frs.^2+(Mpars(2)+Lpars(1))*frs+(Mpars(3)+Lpars(2)))';
right=(Mpars(1)*frs.^2+(Mpars(2)+Rpars(1))*frs+(Mpars(3)+Rpars(2)))';


%other properties needed--------------------------------
spotpos=ThisRep.PosKyTracCom.trackpos_ext;



%define some more arrays for upcoming analysis
lfr=length(frs); 
spotexcess=zeros(lfr,1);
nofits=zeros(lfr,3); 
prefits=zeros(lfr,7);  %parameters of 2x 1D fit per frame
finalfits=zeros(lfr,7); %parameters of 2D fit per frame
spotnos=zeros(lfr,1);   %'number of spots' per frame



%%-----------------------------------------------------------------
% Run a loop on the bacterial frames & positions to fetch
%corresponding images %later to be expand to twoD fit)
dif1=[]; dif2=[];     
for k=1:lfr  %for all frames
lfr-k;
if show==1,close all; end

%get input for 'this' frame --------------------------------
% fitpresets.darklevel=ThisRep.FluoPropsGen.darklevel_ext(k);
% fitpresets.freelabellevel=ThisRep.FluoPropsGen.fluospinelevel_ext(k);
thisfr=round(frs(k));
thisleft=round(max([(left(k)) 1]));    
thisright=round(min([(right(k)) initval.kymolength])); 
thischan_FL=fliplr(squeeze(chanstk_FL(:,:,thisfr)));

%create ROI;%---------------------Do 1D and 2D analysis------------- --------------
if thisright>thisleft+6  %if this is a reasonable ROI to fit

FLspot=thischan_FL(initval.kymohwidth-fitpresets.hwiy+1:initval.kymohwidth+fitpresets.hwiy+1,thisleft:thisright)
FL=thischan_FL(initval.kymohwidth-fitpresets.hwiy+1:initval.kymohwidth+fitpresets.hwiy+1,thisleft:thisright);

fluoprops=Processing_Fluorescence_PatternAnalysis(FLspot); %Get some general pattern properties from this roi

[f1D,f2D]=Processing_stepwisetriplePeakFit(FLspot,fluoprops,fitpresets,plotit); %pre-fit, finalfit

f1D=Processing_EvalGaussFits(f1D,fluoprops,FLspot,initval);
f2D=Processing_EvalGaussFits(f2D,fluoprops,FLspot,initval);

f1D.X0=f1D.X0-1+thisleft;  %transform to absolute coordinates
f1D.X1=f1D.X1-1+thisleft;  %transform to absolute coordinates
f2D.X0=f2D.X0-1+thisleft;  %transform to absolute coordinates
f2D.X1=f2D.X1-1+thisleft;  %transform to absolute coordinates

[I1, I2,I_all]=Get_spotarea_Intensity(f1D,FLspot,fluoprops,initval,thisleft);
%Get Intensities per spot and total (overlap corrected); free-label level and background are subtracted

else   
[f1D,f2D,I1, I2,I_all]=Set_as_bad_point;
end


%collect data points---------
F1X0(k)=f1D.X0;
F1X1(k)=f1D.X1;
F1Y0(k)=f1D.Y0;
F1Y1(k)=f1D.Y1;
F1B(k)=f1D.contentcyto;
F1S1(k)=f1D.contentspot1;
F1S2(k)=f1D.contentspot2;
F1A(k)=f1D.contenttotal;
F1OK1(k)=f1D.Spot1OK;
F1OK2(k)=f1D.Spot2OK;

F2X0(k)=f2D.X0;
F2X1(k)=f2D.X1;
F2Y0(k)=f2D.Y0;
F2Y1(k)=f2D.Y1;
F2B(k)=f2D.contentcyto;
F2S1(k)=f2D.contentspot1;
F2S2(k)=f2D.contentspot2;
F2A(k)=f2D.contenttotal;
F2OK1(k)=f2D.Spot1OK;
F2OK2(k)=f2D.Spot2OK;

nofits(k,:)=[I1,I2,I_all];  %estimated 'no-fit' based on area estimates via f1


%--------------
%'pictures' ; %'pictures';% 'geometry'; 'pictures'; %'geometry'; %'pictures';%'geometry'; %'pictures';% 'geometry'  %'pictures'
if show==1, wht='pictures'; else wht=0;end
switch wht
case 'pictures'
% %plot pictures
% pcolor(thischan_FL); shading flat; colormap hot; title('FL');
% axis equal,
whitebg('white');
subplot(2,1,1); 
lbl=num2str(k);
title(strcat('Replicationcluster',lbl));
subplot(2,1,2); %figure; 
pcolor(FLspot); shading flat; colormap hot; title('FL'); 
%axis([1 2*fitpresets.hwi+1 1 2*fitpresets.hwiy+1]), 
title(strcat('frameno',num2str(k)));hold on
whitebg('white');
plot(f1D.X0+0.5,f1D.Y0+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','b');
plot(f1D.X1+0.5,f1D.Y1+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','g');
plot(f2D.X0+0.5,f2D.Y0+0.5,'o','MarkerFaceColor', 'b', 'MarkerSize', 12, 'MarkerEdgeColor','b');
plot(f2D.X1+0.5,f2D.Y1+0.5,'o','MarkerFaceColor', 'g', 'MarkerSize', 12, 'MarkerEdgeColor','g'); 
%if spotno==2 
if (~isnan(f2D.Spot1OK) & ~isnan(f2D.Spot2OK))
plot([f2D.X0 f2D.X1]+0.5, [f2D.Y0 f2D.Y1]+0.5, '-');
end
hold off;
%[~]=ginput(1);
pause(1.0);

if mov
%hold on 
%outname=strcat(filename, int2str(j));
%saveas(gcf,outname,'jpg')
%F(k) = getframe;
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

%finally, store in cluster
areasums.I1=nofits(:,1);
areasums.I2=nofits(:,2);
areasums.Iall=nofits(:,3);


prefits.X0=F1X0;
prefits.X1=F1X1;
prefits.Y0=F1Y0;
prefits.Y1=F1Y1;
prefits.contentcyto=F1B;
prefits.contentspot1=F1S1;
prefits.contentspot2=F1S2;
prefits.contenttotal=F1A;
prefits.spot1OK=F1OK1;
prefits.spot2OK=F1OK2;

finalfits.X0=F2X0;
finalfits.X1=F2X1;
finalfits.Y0=F2Y0;
finalfits.Y1=F2Y1;
finalfits.contentcyto=F2B;
finalfits.contentspot1=F2S1;
finalfits.contentspot2=F2S2;
finalfits.contenttotal=F2A;
finalfits.spot1OK=F2OK1;
finalfits.spot2OK=F2OK2;

end

%Test plot menu
if 1
    close(gcf);
    scrsz = get(0,'ScreenSize');
    figure('Position',[100 100 scrsz(3)/1.3 scrsz(4)/1.3]);
    subplot(2,2,1);
    plot(frs,left, 'k-'); hold on;
    plot(frs,right, 'k-');
    plot(frs, spotpos, 'k-o');
    plot(frs,prefits.X0.*prefits.spot1OK, 'b-o');
    plot(frs,prefits.X1.*prefits.spot2OK, 'r-o');
    %[~]=ginput(1); 
    pause(0.1);
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


function [I1, I2,I_all, N1,N2,Nall]=Get_spotarea_Intensity(f1,FLspot,fluoprops,initval,thisleft);
%This function obtains the summed intensity around predetermined maxima
%points; %prefits: x0,x1,y0,y1, b0,N0, N1
%first peak is brightest


% fluo_roi = 
%         darklevel: 162 background level
%      (all values below are background corrected using 'darklevel')
%     signalcontent: 314049 sum(pixel values-background)
%     sumspinelevel: 6742 median of y-summed pixels
%     maxspinelevel: 1355 median of y-max pixels (representative of
%     cytoplasm level)
%      spotscontent: 219751 summed values of pixels above sumspinelevel
%           peakval: 15612  maximum main peak
%             peaky: 10     position peak
%             peakx: 12
%FLspot: image with spots between edges of bacterium


hw=initval.sumHW; %area width is (hw*2+1);

if 0 %plot menu--------------------------------------------------
close all
pcolor(FLspot); shading flat;

pcolor(FLspot); shading flat; colormap hot; title('FL');hold on
whitebg('white');
plot(f1.X0-thisleft+1+0.5,f1.Y0+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','b');
plot(f1.X1-thisleft+1+0.5,f1.Y1+0.5,'o', 'MarkerSize', 12, 'MarkerEdgeColor','g');
[~]=ginput(1);
end

[r,c]=size(FLspot)
%border patrol------------------
x1=max(f1.X0-thisleft+1, 1+hw); x1=min(x1, c-hw);
y1=max(f1.Y0, 1+hw); y1=min(y1, r-hw);
x2=max(f1.X1-thisleft+1, 1+hw); x2=min(x2, c-hw);
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

%summed intensities; correcting for bacterium width
I1=sum(area1(:))-fluoprops.darklevel*N1-2*hw/7*fluoprops.sumspinelevel;  %assuming bacterium has width 5; spot 2.5
I2=sum(area2(:))-fluoprops.darklevel*N2-2*hw/7*fluoprops.sumspinelevel;
I3=sum(area3(:))-fluoprops.darklevel*N3-2*hw/7*fluoprops.sumspinelevel;

if ~isempty(I3)
I_all=I1+I2-I3;
else
I_all=I1+I2;
end
end

function [f1D,f2D,I1, I2,I_all]=Set_as_bad_point;
%Dummy function to set NaN vals
f1D.X0=NaN;
f1D.X1=NaN;
f1D.Y0=NaN;
f1D.Y1=NaN;
f1D.contentcyto=NaN;
f1D.contentspot1=NaN;
f1D.contentspot2=NaN;
f1D.contenttotal=NaN;
f1D.Spot1OK= NaN;
f1D.Spot2OK= NaN;

f2D.X0=NaN;
f2D.X1=NaN;
f2D.Y0=NaN;
f2D.Y1=NaN;
f2D.contentcyto=NaN;
f2D.contentspot1=NaN;
f2D.contentspot2=NaN;
f2D.contenttotal=NaN;
f2D.Spot1OK= NaN;
f2D.Spot2OK= NaN;
f2D.Spot1OK= NaN;
f2D.Spot2OK= NaN;

I1=NaN;
I2=NaN;
I_all=NaN;

fluo.darklevel=NaN;
fluo.signalcontent=NaN;
fluo.maxspinelevel=NaN;
fluo.sumspinelevel=NaN;
fluo.spotscontent=NaN;
fluo.peakval=NaN;
fluo.peaky=NaN;
fluo.peakx=NaN;
end

