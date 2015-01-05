
function pc=GetParticles(im,initval);
    %This function finds significant spots in an image. We assume the image
    %is pre-processed to units of sigma, where sigma is the local noise
    %level 
    edz=initval.edge;
    im = smooth(im,1);                               %smooth a bit to bring in implicit integration step (and limit number of maxima)
    [r,c]=size(im);    
    [mc,max_vals] = findmaxima(im(edz:r-1-edz,edz:c-1-edz));  %stay away from the edges
    mc = mc + edz;                            %correct for shift caused by edge truncation in previous line
    mc = round(mc);                       %make coords integer 
    flag=Outlier_Flag(max_vals,initval.particletreshold,0.8,'positive',0);
    ind=find(flag==0);                      %the outliers are particles!    
    pc = [mc(ind,1) mc(ind,2) max_vals(ind) 0*max_vals(ind)]; %last column for later addition frame no
    dum=1;
end