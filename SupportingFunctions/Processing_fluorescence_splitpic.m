function [FL_spot,FL_cyto,FL_residu]=Processing_fluorescence_splitpic(FL,fluo);
%This function returns 'spot-only' and 'cytoplasm only' pictures, based on simple
%assumption of the fluorescent pattern of the bacterium. these are intended
%to be used for fitting of the spots
%-------------------------------------
% fluo = 
        % area_bac: 118
        % area_spot: 23
        % content_cytoplasm1: 196929.939867696
        % content_cytoplasm2: 225558.998548967
        % content_signal: 299074.998548967
        % content_spots1: 102145.058681271
        % content_spots2: 73516
        % content_total: 760405
        % curve_medianofmax: [1x17 double]
        % curve_medianofsum: [1x17 double]
        % level_dark: 749.264705882353
        % level_fluotreshold: 906.764710624292
        % level_medianofmax: 2837.23528937571
        % level_medianofmax_yposcurve: [1x17 double]
        % level_medianofsum: 13027.5293406356
        % level_peak: 13839
        % noise_dark: 78.7500023709693
        % peak_xpos: [1x17 double]
        % peak_ypos: [1x17 double]
        % ratio_FS: 0.245811252550966
        % ratio_SN: 27.855834798822
        % wherebac: [118x1 double]
        % wheredark: [238x1 double]
        % wherefluo: [306x1 double]
        % wherespot: [23x1 double]
%Note: all 'levels' beyond dark are calculated from this dark level
%------------------------------------------


if nargin<2 %TEST MODUS; using Charls Database of single-focus-images
    close all; clear all;
    pth='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\2013_08_14 FociEval\ImageDatabase\SingleFocus\';
    %pth='D:\jkerssemakers\My Documents\BN_ND_ActiveProjects\BN_ND11_CharlBacterialReplication\2013_08_14 FociEval\ImageDatabase\NoFoci\';
    imagenames=dir(strcat(pth ,'*.tif'));
    nm=strcat(pth,imagenames(5).name);
    im0=double(imread(nm));  %load an image; including borders         
    [r,c]=size(im0);        %Crop to nonzero area
    [X,Y]=meshgrid(1:c,1:r);
    sel=find(im0~=0);
    lor=min(Y(sel));             hir=max(Y(sel)); 
    loc=min(X(sel));             hic=max(X(sel));
    FL=im0(lor:hir,loc:hic); 
    [fluo,modelpic]=Processing_Fluorescence_PatternAnalysis(FL);
end

    %first, we define a backbone of this bacterium:

    [r,c]=size(FL);
    bb0_x=[1:c]; axy=[1:c]; [X,Y]=meshgrid(bb0_x,axy);
    bb0_y=fluo.level_medianofmax_yposcurve;
    bb0_I=fluo.curve_medianofmax;

    %select the higher part; build a picture with one backbone contour that
    %with an intensity of the median (y-summed) fluorescence counts; these
    %counts represent the cytoplasmic counts
    sel=find(bb0_I>fluo.level_fluotreshold); 
    lox=min(sel); hix=max(sel);
    bb1_x=bb0_x(lox:hix); 
    bb1_y=bb0_y(lox:hix); 
    bb1_I=bb0_I(lox:hix);
    FL_bb=0*FL;
    for i=1:length(bb1_x)
        FL_bb(bb1_y(i),bb1_x(i))=fluo.level_medianofsum;
    end

    %Blur this backbone picture to form a 'cytoplasm' picture  
    hpsf=1.7;
    squ=5;
    [xk,yk]=meshgrid(1:squ,1:squ);
    radius=((xk-squ/2).^2+(yk-squ/2).^2).^0.5; %'radial picture'
    blurkernel=1/(2*pi*(hpsf.^2))*exp(-radius.^2/(2*hpsf.^2));
    FL_cytobar=imfilter(FL_bb,blurkernel);

    %use this cytoplasmic picture to render a 'spot only' picture in which only
    %spots are left above a flat noise background
    FL_spotn=FL-FL_cytobar;
    FL_spot=0*FL_spotn;    %noise plus spot
    sel=find(FL_spotn>2*fluo.level_fluotreshold);
    FL_spot(sel)=FL_spotn(sel)-fluo.level_fluotreshold; %only spot

    FL_cyto=0*FL;
    sel=find(FL>fluo.level_fluotreshold);
    FL_cyto(sel)=FL(sel)-fluo.level_fluotreshold;
    FL_cyto=FL_cyto-FL_spot;                           %only cyto
    FL_residu=FL-FL_cyto-FL_spot;


if nargin<1
    close all;
    figure;
    %plot trickses
    scalebar=max(FL(:));
    subplot(2,4,1);
    FL(1,1)=scalebar;
    pcolor(FL); shading flat; hold on;
    plot(bb1_x,bb1_y,'w-', 'MarkerFaceColor', 'w'); hold on

    subplot(2,4,2);
    FL_bb(1,1)=scalebar-fluo.level_dark;
    pcolor(FL_bb); shading flat; hold on;

    subplot(2,4,3);
    FL_cytobar(1,1)=scalebar-fluo.level_dark;
    pcolor(FL_cytobar); shading flat; hold on;

    subplot(2,4,4);
    FL_spotn(1,1)=scalebar;
    pcolor(FL_spotn); shading flat; hold on;

    subplot(2,4,5);
    FL_spot(1,1)=scalebar-fluo.level_dark;
    pcolor(FL_spot); shading flat; hold on;

    subplot(2,4,6);
    FL_cyto(1,1)=scalebar-fluo.level_dark;
    pcolor(FL_cyto); shading flat; hold on;
    
    subplot(2,4,7);
    FL_residu(1,1)=scalebar
    FL_cyto(1,1)=scalebar-fluo.level_dark;
    pcolor(FL_residu); shading flat; hold on;
end
