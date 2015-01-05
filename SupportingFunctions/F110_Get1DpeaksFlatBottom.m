function mxi=F110_Get1DpeaksFlatBottom(data,sigs);
% peakfinding routine
%Find peaks in a profile 'data'. A peak is a local maximum that is higher
%than a certain treshold.  This treshold is found assuming a 'flat' bottom, i.e. peaks sticking out of a reasonably flat base level.
%This base level is used to find a proper measure for the standard deviation 

%input: 1D data containing enough points to perform statistics; 
%output: indices of the accepted peaks
%Jacob Kers 2013
%--------------------------------------------------------------------------
plotit=0;


if nargin<1  %For testing puroposes
close all
data=Demodata;
sigs=2;
end

md=median(data);
sigma1=std(data);

[flag,cleandata]=F_10000_Outlier_Flag(data,3,0.7,'positive',0);
sigma2=std(cleandata);
md2=median(cleandata);

localmaxidx=find(data(2:end-1)>data(1:end-2) &data(2:end-1)>data(3:end))+1;  %local maxes 1D
sel=find(data(localmaxidx)>md2+sigs*sigma2);
mxi=localmaxidx(sel);

if plotit;
    plot(data,'o-'); hold on;  %For testing puroposes
    stem(mxi,data(mxi), 'ko', 'MarkerFace', 'k'); %vertical lines at peak positions
end


function data=Demodata
data=[   609
         860
         894
         898
         587
         617
         807
         629
         577
         531
        1362
        1152
        1135
        1738
        4247
        9516
       10652
        6829
        3294
        1523
         977
        2158
         910
         502
         840
        1068
         530
         855
         810
         787
        1124
        1296
        3502
        6612
        9071
        8176
        6293
        4184
        1589
        1724
        1718
        1447
         736
         781
         948
         603
         429
         587
         630
         581
         748
        1040
        1151
         642
        1086
         876
        1907
        4324
        5548
        6736
        4446
        3718
        2919
        2089
        1874
        1963
        1621
        1144
        1417
         928
        1024
        1634
         834
        1907
        1729
        3508
        9138
       13053
        6597
        2018
        1303
        1214
         985
         793
        1140
        1694
        1376
        1682
        1685
        3641
        4964
        7704
        8733
        6773
        9497
        4653
        2803
        1904
        2006
        2023
        1343
        1730
        1553
        1098
        1892
        1635
        2089
        2849
        4702
        8142
        7245
        5983
        2935
        1727
        1581
         730
        1405
         861
         838
         576
         990
         381
         453
         532
         469
         587
         603
         217
        1262
         324
         386
         578
         420
         437
         371
         495
         583
         402
         185
          16
           0
           0
           0
           0
           0
           0
           0
           ];

