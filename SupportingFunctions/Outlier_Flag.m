
function [flag,cleandata]=Outlier_Flag(data,tolerance,sigchange,how,sho,binz);
%this function calculates the standard deviation and average of a chosen column of the data; Then it
%throws out the rows that contain in that column a value considered
%unreasonable. This is repeated until the new sigma does not change much
%anymore
%output: positions of outliers
    %figure;
    sigma=1E20;            %at start, use a total-upper-limit 
    ratio=0;
    ld=length(data);
    flag=ones(ld,1);  %at start, all points are selected
    cleandata=data;
    while ratio<sigchange     %if not too much changes anymore; the higher this number the less outliers are peeled off.
        sigma_old=sigma;
        selc=find(flag==1);
        data(selc); 
        ls=length(selc);
        av=nanmedian(data(selc));       %since we expect skewed distribution, we use the median iso the mea     
        sigma=nanstd(data(selc));
        ratio=sigma/sigma_old;
        switch how
            case 'positive',  flag=(data-av)<tolerance*sigma;     %adjust outlier flags
            case 'all',  flag=abs(data-av)<tolerance*sigma;     %adjust outlier flags  
        end
        
    end
    cleandata=data(selc); 
    hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
    sthst=hist(cleandata,hx);
    
    if sho==1
            figure;
            bar(hx,sthst);
            title('Histogram');
            dum=ginput(1);
            pause(0.5);     
        end
end