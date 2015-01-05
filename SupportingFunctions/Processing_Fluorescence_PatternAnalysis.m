function [fluo,modelpic]=Processing_Fluorescence_PatternAnalysis(FL);
%This function returns some fluorescent properties, based on simple
%assumption of the fluorescent pattern of the bacterium

%-------------------------------------
% fluo = 
% 
%               area_bac: 118
%              area_spot: 23
%     content_cytoplasm1: 196929.939867696
%         content_signal: 299074.998548967
%         content_spots1: 102145.058681271
%          content_total: 760405
%      curve_medianofmax: [1x17 double]
%      curve_medianofsum: [1x17 double]
%     level_cytotreshold: 906.764710624292
%             level_dark: 749.264705882353
%      level_medianofmax: 2837.23528937571
%      level_medianofsum: 13027.5293406356
%             level_peak: 13839
%             noise_dark: 78.7500023709693
%              peak_xpos: [1x17 double]
%              peak_ypos: [1x17 double]
%               ratio_FS: 0.245811252550966
%               ratio_SN: 27.855834798822
%               wherebac: [118x1 double]
%              wherefluo: [306x1 double]
%              wherespot: [23x1 double]
%Note: all 'levels' beyond dark are calculated from this dark level
%------------------------------------------


if nargin<1 %TEST MODUS; using Charls Database of single-focus-images
    close all
    pth='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\2013_08_14 FociEval\ImageDatabase\SingleFocus\';
    %pth='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\2013_08_14 FociEval\ImageDatabase\NoFoci\';
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
fluo.level_fluotreshold=fluotreshold;
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
fluo.content_spots=sum(fluosumcurve(sel)-fluo.level_medianofsum);

%b) As a second approach to counting spot content, we use the 'maxima' median
%level and the full picture. This 2D-approach may give us a sharper distinction between
%'spot' and 'no-spot' area
[fluomaxcurve,maxidx]=max(FLbc);


fluo.level_medianofmax=median(fluomaxcurve);  %level represents backbone level bacterium
fluo.level_medianofmax_yposcurve=maxidx;  %line showing position of maxima, spanning whole picture
fluo.curve_medianofmax=fluomaxcurve;
fluo.wherespot=find(FLbc>fluo.level_medianofmax);
spotarea=length(fluo.wherespot);  %amount of 'spot pixels'.
        %'shave off' above-median-counts
fluo.area_spot=spotarea;


%as a convenient number to use, we also measure the bacterium area (taking
%a stricter selecection of pixels than for the counts
edgetreshold=fluo.level_dark+8*fluo.noise_dark;
fluo.wherebac=find(FL>edgetreshold);
fluo.area_bac=length(fluo.wherebac);

%make a picture representing cytosol, and spot area, to show if the selected
%area cover realistic parts; store these indices for later use
fluo.wherefluo=find(FL>fluotreshold);
fluo.wheredark=find(FL<=fluotreshold);
modelpic=0*FL;
modelpic(fluo.wherefluo)=1;
modelpic(fluo.wherebac)=2;
modelpic(fluo.wherespot)=3;

%find global max as crude spot location
fluo.level_peak=max(FL(:));
[~,fluo.peak_xpos]=max(max(FL));
[~,fluo.peak_ypos]=max(max(FL'));


fluo.content_cytoplasm=fluo.content_signal-fluo.content_spots;

%Last, detrmine some handy ratios
fluo.ratio_FS=fluo.content_spots/fluo.content_signal;
fluo.ratio_SN=mean(fluomaxcurve)/(2*fluo.noise_dark);

else 
fluo.area_bac=0;
fluo.area_spot=0;
fluo.content_cytoplasm=0;
fluo.content_signal=0;
fluo.content_spots=0;
fluo.content_total=0;
fluo.curve_medianofmax=1;
fluo.curve_medianofsum=1;
fluo.level_fluotreshold=0;
fluo.level_dark=0;
fluo.level_medianofmax=0;
fluo.level_medianofmax_yposcurve=0;
fluo.level_medianofsum=0;
fluo.noise_dark=0;
fluo.ratio_FS=0;
fluo.ratio_SN=0;
fluo.level_peak=0;
fluo.peak_xpos=0;
fluo.peak_ypos=0;
fluo.wherebac=[];
fluo.wheredark=[];
fluo.wherefluo=[];
fluo.wherespot=[]; 
end

fluo=orderfields(fluo);


if nargin<1
fluo
subplot(2,2,1); pcolor(FL); shading flat; 
title('original image');
subplot(2,2,2); pcolor(modelpic);
title('regions: background-fluorescence-edge-spots');
subplot(2,2,3); 
    plot(fluosumcurve, 'k-o');hold on; 
    dummy=0*fluosumcurve;
    plot(dummy+fluo. level_medianofsum,'r-');
    axis([0 length(sum(FL)) 0 1.2*max(fluosumcurve)]);
    title('sums of X-sections');
    xlabel('x-position (pixel units)');
    ylabel('Summed Intensity')
subplot(2,2,4); 
    plot(fluomaxcurve, 'k-o');hold on; 
    dummy=0*fluosumcurve;
    plot(dummy+fluo. level_medianofmax,'r-');
    axis([0 length(sum(FL)) 0 1.2*max(fluomaxcurve)]);
    title('maxima of X-sections');
    xlabel('x-position (pixel units)');
    ylabel('Max.Intensity')
[~]=ginput(1);
close(gcf);
end

end %if test