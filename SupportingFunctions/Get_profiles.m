
function bleachcurves=Get_profiles(aa,pc, initval);

%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%III Get 'profiles' : vertical slices through the stack of images per particle
%location; each profile reflects time-dependent intensity changes of that
%particle. Each intensity point is taken from a
%small (2x2-4x4) area around the particle location.
 %Get Z-profiles by finding local maxes around xy positions of detected particles   
    [r,c,d]=size(aa);        plns=1:d;  lm=length(pc); 
    all_plns=zeros(d*lm,2); %concatenated bleach curves
    av_all_plns=zeros(d*lm,2); %concatenated bleach curves
    
    all_plns_cols=zeros(lm,d); %column-grouped bleach curves
    figure;
    ab=dip_array(aa);
    for i=1:lm 
        
           
            pos=[pc(i,1), pc(i,2)];
            %profile=Get_averagedZ_profile_useLocalInfo(aa,pos, initval);
            [profile,avprofile]=Get_Zprofile2DGaussianFit(ab,pos,initval);
            profile=profile(1:d); 
            avprofile=avprofile(1:d);
            
            
            all_plns((i-1)*d+1:(i-1)*d+d,1)=[1:d];      %concatenated bleach curves
            all_plns((i-1)*d+1:(i-1)*d+d,2)=profile;    %column-grouped bleach curves
            
            av_all_plns((i-1)*d+1:(i-1)*d+d,1)=[1:d];      %concatenated bleach curves
            av_all_plns((i-1)*d+1:(i-1)*d+d,2)=avprofile;    %column-grouped bleach curves
            
            
            all_plns_cols(i,:)=profile;
           
            %plot(plns,profile+0*(i-1)*1000, '-');  hold on; % axis([0 d 0 2000]);
            xs=i*d+d-1; ra=max(av_all_plns(1:xs,2));
            plot(all_plns(1:xs,2)); hold on  
            plot(av_all_plns(:,2),'r'); 
            axis([0 xs -2000 ra])
%             close all;
%             subplot(2,1,1);
%             plot(plns,profile,'-'); hold on
%             plot(plns,smooth(profile,5),'r');
%            
%             axis([0 d 0 ra]);
%             binz=50;
%             hx=(0:ra/binz:ra);   %make an axis
%             sthst=hist(all_plns,hx);           %make a histogram
%             subplot(2,1,2);
%             bar(hx,log10(sthst));
%             title('Histogram');
             pause(0.1); 
    end
    bleachcurves=all_plns_cols;
    avcurve=mean(all_plns_cols);
    figure;
    plot(avcurve);
    dum=1;
    lbl=(strcat(initval.datapath,'bleachcurves.txt'));
    dlmwrite(lbl,all_plns(:,:),'delimiter',' ', 'precision', 6);
end

