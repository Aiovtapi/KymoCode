function pc=Find_Particles_JKCM(im,initval);
%This code is intended to collect bleaching curves of particle image 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%NOTES: 
%'dipimage' should be installed and started ('dipstart') to make this
%code work
%TIP: run timesseries reader once, then save the
%workspace including 'a3' ;then comment out this section and run 
%'dum=AD2011_11_projectBR_BleachcurveAnalyzer(a3)' from the command line; 
%saves lots of %time with repaeated runs
%Jacob Kerssemakers 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initval.radialselect=120;     % pixel radius wherin particles  will be used for density calculation
initval.standardarea=20;       %particle counts per ..X.. microns
initval.nmperpixel=160;

%initval.resultpath='D:\jkerssemakers\My Documents\BN_Active_Projects\BN_ND11_BacterialReplicationCharl\';
%initval.datastring='D:\jkerssemakers\My Documents\MATLAB\2011_BR BeadCountPack\ParticleCountTest\bleaching_test_1971t*.tif'
initval.backsize=25;       %number of regions used for creating initial background image
initval.minproximity=5;    %particles should be more apart than this value
%initval.minproximity=25;    %helps for bright beads
initval.square=6;
%initval.sigma=2.5;          %treshold for detecting particles
warning off all
%close all; 

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(b)Make estimate of 'mean field' image  by low-pass
%filtering the maximum-projected image 
    sq=20;
    %sq=50;  %works better for bright beads
    h = ones(sq,sq) /sq^2; 
    mean_im = imfilter(im,h);                     
    %dipshow(4,mean_im, [0 2]);
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(b)Make estimate of 'local noise' image ; 
%Assume this noise to stem from shot noise; then illumination would be
%square of this:     
    flat_im=(im)-mean_im;                          %substract 'mean field'    
    [noise_im,illum_im]=GetLocalNoiseImage(initval,flat_im) ;
    %dipshow(5,noise_im, [0 200]); 
    %dipshow(6,illum_im, [0 1.5]);
   
 %finally, make an image to be used for particle detection
 particle_im=flat_im./noise_im;
 %dipshow(8,particle_im, [0 5]);
 
 
 %III. Using these images, we find particles by searching what sticks out of the noise   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get particles; This is done by detecting local maxima in the max-projected image
%sufficiently above the background; particles too close together are
%rejected.   
     pc=GetParticles(particle_im,initval);  %Get peaks that are sufficiently above the background
     pc=Do_proximity_test(pc,initval);               % Reject peaks too close to a brighter one
     result=strcat('detected:  ',num2str(length(pc)) ,' particles')
 
     
     
     
end


function back_sm=GetBackgroundImage(initval,b)
    %This function generates a background images by first taking the medians of
    %sub_regions and next smoothing the result
    back_image = newim(b);              %create empty background image
    
    [r,c]=size(back_image);
    for i = 0:initval.backsize:r-initval.backsize
        for j = 0:initval.backsize:c-initval.backsize
            back_image(i:i+initval.backsize-1,j:j+initval.backsize-1) = median(b(i:i+initval.backsize-1,j:j+initval.backsize-1));
        end
    end
    back_sm = smooth(back_image,1.5*initval.backsize);
    %dipshow(1,back_sm,[550,2500]);
end

function [noise_image,illum_image]=GetLocalNoiseImage(initval,b)
    %This function produces a 'noise image' , where every pixel value
    %represents a local estimate of 1 sigma of a presumed gaussian intensity noise
    %distribution (it is poissonian, however)
     %the image is divided in 12x12 areas; 
    %this assumes a reasonable smooth varying illumination (and therefore, noise level)
    %'reasonable' would be a gaussian with sigma=~imagesize/2 or larger   
     noise_image = 0*b;              %create empty background image
    [r,c]=size(noise_image);
    sba=ceil(r/16);                  
    for i = 1:sba:r-sba+1
        for j = 1:sba:c-sba+1
            sub_area=b(i:i+sba-1,j:j+sba-1);
            noise=Get_Noise_from_pairedneighbours(sub_area);
            noise_image(i:i+sba-1,j:j+sba-1) = noise;
        end
    end
    noise_image = smooth(noise_image,1.5*sba);
    
    bf=noise_image.^2;
    illum_image=bf/max(max(bf));       %normalized illumination image (assuming poissonian noise)
end


function noise=Get_Noise_from_pairedneighbours(arr);
%This function estimates the noise in an image by pair-wise difference
%(thus excluding slow gradients like background
    [r,c]=size(arr);
    difim=arr(:,2:c)-arr(:,1:c-1);
    noise=mean(std(difim))/2^0.5;
end
     

function pc=GetParticles(im,initval);
    %This function finds significant spots in an image. We assume the image
    %is pre-processed to units of sigma, where sigma is the local noise
    %level
    
    im = smooth(im,1);                               %smooth a bit to bring in implicit integration step (and limit number of maxima)
    %dipshow(11,im,[0,5]);
    [r,c]=size(im);    
    [mc,max_vals] = findmaxima(im(10:r-11,10:c-11));  %stay away from the edges
    mc = mc + 10;                            %correct for shift caused by edge truncation in previous line
    mc = round(mc);                       %make coords integer 
    flag=Outlier_Flag(max_vals,initval.sigma,0.9);
    ind=find(flag==0);                      %the outliers are particles!    
    pc = [mc(ind,1) mc(ind,2) max_vals(ind)];    
    dum=1;
end


function pc2=Do_proximity_test(pc,initval)
%This function calculates for each position if a too nearby maximum is
%larger or not. If so, it is discarded. JacobKers11
    pc2=0*pc;
    [lp,~]=size(pc);  c=0;
    for i=1:lp
        x0=pc(i,1); y0=pc(i,2); I=pc(i,3);    %xyI coordinates
        dist=((pc(:,1)-x0).^2+(pc(:,2)-y0).^2).^0.5;
        sel=find((dist<initval.minproximity)&dist~=0);
        if ~isempty(sel);       %close neighbours
            Iprox=pc(sel,3); 
            if I>2*max(Iprox)     %is it -clearly- the local winner?
                c=c+1;
                pc2(c,:)=pc(i,:);
            end
        else
          c=c+1;%isolated maximum
          pc2(c,:)=pc(i,:);    
        end
    end
    pc2=pc2(1:c,:);
end

     
 
function flag=Outlier_Flag(data,tolerance,sigchange);
%this function calculates the standard deviation and average of a chosen column of the data; Then it
%throws out the rows that contain in that column a value considered
%unreasonable. This is repeated until the new sigma does not change much
%anymore
%output: positions of outliers
    %figure;
    binz=25;
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
        flag=(data-av)<tolerance*sigma;     %adjust outlier flags
        hx=(min(data(selc)):(range(data(selc)))/binz:max(data(selc)));   %make an axis
        sthst=hist(data(selc),hx);       
%         figure;
%         bar(hx,sthst);
%         title('Histogram');
%         dum=ginput(1);
%         pause(0.2);
%         close(gcf);
    end
    cleandata=data(selc); 
    hx=(min(cleandata):(range(cleandata))/binz:max(cleandata));   %make an axis
    sthst=hist(cleandata,hx);
end

 
   

    