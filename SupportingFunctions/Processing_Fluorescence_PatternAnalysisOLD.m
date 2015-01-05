function [fluo,modelpic]=Processing_Fluorescence_PatternAnalysis(FL);
%This function returns some fluorescent properties, based on simple
%assumption of the fluorescent pattern of the bacterium
% fluo = 
% 
%               area_bac: 157
%              area_spot: 24
%     content_cytoplasm1: 356080.232526967
%         content_signal: 416051.793024526
%         content_spots1: 59971.5604975588
%          content_total: 1426534
%      curve_medianofmax: [1x31 double]
%      curve_medianofsum: [1x31 double]
%     level_cytotreshold: 1157.35514397638
%             level_dark: 931.129032258065
%           level_maxval: 6394.64485602362
%      level_medianofmax: 3726.64485602362
%      level_medianofsum: 13311.8037042126
%             noise_dark: 113.113055859158
%               ratio_FS: 0.144144458702101
%               ratio_SN: 15.6263856637758
%     
%FLspot: image with spots between edges of bacterium

if nargin<1 %TEST MODUS; using Charls Database of single-focus-images
    close all
    %pth='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\2013_08_14 FociEval\ImageDatabase\SingleFocus\';
    pth='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\2013_08_14 FociEval\ImageDatabase\NoFoci\';
    imagenames=dir(strcat(pth ,'*.tif'));
    nm=strcat(pth,imagenames(1).name);
    im0=double(imread(nm));  %load an image; including borders         
    [r,c]=size(im0);        %Crop to nonzero area
    [X,Y]=meshgrid(1:c,1:r);
    sel=find(im0~=0);
    lor=min(Y(sel));             hir=max(Y(sel)); 
    loc=min(X(sel));             hic=max(X(sel));
    FL=im0(lor:hir,loc:hic); 
end


[r,c]=size(FL);
if c>1
fluo.content_total=sum(FL(:));


%1 First, we make estimates on intensity and noise of the background (far outside
%the bacterium).'local background' is defined as the average of the outer two image lines'
fluo.level_dark=(mean(mean(FL(1:2,:)))+mean(mean(FL(r-1:r,:))))/2;


%2) estimate dark (camera) noise via the spread in the difference between
%neighbouring pixels
diftop=FL(1:5,2:end)-FL(1:5,1:end-1); 
difbot=FL(r-4:r,2:end)-FL(r-4:r,1:end-1);
dif=[diftop(:);  difbot(:)];
fluo.noise_dark=std(dif)/2^0.5;

%define as 'fluorescence' those pixels sufficiently above the darklevel.
%Note this may not be representative for the outline of the bacterium,
%since there is some blurring and we want to measure all fluorescence
fluotreshold=fluo.level_dark+2*fluo.noise_dark;
FLbc=FL-fluotreshold;
backsel=find(FL<fluotreshold);
fluo.level_cytotreshold=fluotreshold;
FLbc(backsel)=0;



%Use these pixels to make two 1D-curves along the x direction (length of bacterium): 
%-one that contains all the counts; 
%one that contains the maxima along the y-direction.
%these curves are used to find 'median' levels, representative for the
%cytoplasmic signal level
%a) First, the summed intensities (minus background)
fluosumcurve=sum(FLbc);
fluo.level_medianofsum=median(fluosumcurve);
fluo.curve_medianofsum=fluosumcurve;
fluo.content_signal=sum(fluosumcurve);       %total represents total counts
sel=find(fluosumcurve>fluo.level_medianofsum);
fluo.content_spots1=sum(fluosumcurve(sel)-fluo.level_medianofsum);

%b) As a second approach to counting spot content, we use the 'maxima' median
%level and the full picture. This 2D-approach may give us a sharper distinction between
%'spot' and 'no-spot' area

fluomaxcurve=max(FLbc);
fluo.level_maxval=max(fluomaxcurve);
fluo.level_medianofmax=median(fluomaxcurve);  %level represents backbone level bacterium
fluo.curve_medianofmax=fluomaxcurve;
wherespot=find(FLbc>fluo.level_medianofmax);
spotarea=length(wherespot);  %amount of 'spot pixels'.
        %'shave off' above-median-counts
%fluo.content_spots2=sum(FLbc(wherespot))-spotarea*fluo.level_medianofmax;
fluo.area_spot=spotarea;


%as a convenient number to use, we also measure the bacterium area (taking
%a stricter selecection of pixels than for the counts
edgetreshold=fluo.level_dark+8*fluo.noise_dark;
bacareasel=find(FL>edgetreshold);
fluo.area_bac=length(bacareasel);

%make a picture representing cytosol, and spot area, to show if the selected
%area cover realistic parts
fluosel=find(FL>fluotreshold);
modelpic=0*FL;
modelpic(fluosel)=1;
modelpic(bacareasel)=2;
modelpic(wherespot)=3;


fluo.content_cytoplasm1=fluo.content_signal-fluo.content_spots1;

%Last, detrmine some handy ratios
fluo.ratio_FS=fluo.content_spots1/fluo.content_signal;
fluo.ratio_SN=mean(fluomaxcurve)/(2*fluo.noise_dark);

if 0
fluo
fluo.content_spots1/fluo.content_signal
subplot(2,1,1); pcolor(FL); pause(0.1);
subplot(2,1,2); plot(fluomaxcurve,'o-'); axis tight; hold on;
subplot(2,1,2); plot(0*fluomaxcurve+fluo.level_medianofmax,'r-'); axis tight; %hold off;
[~]=ginput(1);
end


else
fluo.area_bac=1;
fluo.area_spot=1;
fluo.content_cytoplasm1=1;
fluo.content_signal=1;
fluo.content_spots1=0;
fluo.content_total=1;
fluo.curve_medianofmax=1;
fluo.curve_medianofsum=1;
fluo.level_cytotreshold=0;
fluo.level_dark=0;
fluo.level_medianofmax=0;
fluo.level_medianofsum=0;
fluo.level_maxval=0;
fluo.noise_dark=0;
fluo.ratio_FS=0;
fluo.ratio_SN=0;
end

fluo=orderfields(fluo);


if nargin<1
fluo
plot(fluo.curve_medianofsum, 'k-o');hold on; 
dummy=0*fluo.curve_medianofsum;
plot(dummy+fluo. level_medianofsum,'r-');
axis([0 length(sum(FL)) 0 1.2*max(fluo.curve_medianofsum)]);
[~]=ginput(1);
close(gcf);
end

end %if test