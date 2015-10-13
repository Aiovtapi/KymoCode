function [flag,cleandata]=F_10000_Outlier_Flag(data,tolerance,sigchange,how,sho);
%this function is meant to find a representative value for a standard
%deviation in a heavily skewed distribution (typically, flat data with
% %peaks). It calculates the standard deviation and average the data;
% Based on these, outliers are determined and excluded for a new calculation
% of average and SD; this is repeated until sigma does not change anymore too much
% . This is repeated until the new sigma does not change much
% %anymore
%output: positions of outliers

%Jacob Kers 2013 and before---------------------------------------------

if nargin<5  %For testing/demo purposes
    close all
    data=Demodata;
    tolerance=2;
    sigchange=0.7;
    how='positive';
    sho=1;
    plot(data,'o-');
end
binz=50;
sigma=1E20;            %at start, use a total-upper-limit 
ratio=0;
ld=length(data);
flag=ones(ld,1);  %at start, all points are selected
cleandata=data;
while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
    sigma_old=sigma;
    selc=find(flag==1);
    data(flag==1); 
    ls=length(selc);
    av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
    sigma=nanstd(data(selc));
    ratio=sigma/sigma_old;
    switch how
        case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
        case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
    end
    %plot menu------------------  
    if sho==1
        cleandata=data(selc); 
        hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
        sthst=hist(cleandata,hx);
        bar(hx,sthst);
        title('Histogram');
        dum=ginput(1);
        pause(0.5);     
    end
    %---------------------------- 
end
cleandata=data(selc); 
hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
sthst=hist(cleandata,hx);



function data=Demodata;
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

