
clc;
%restoredefaultpath;
%clc;
clear all;
% close all;
% clear functions;
% bdclose('all');
% fclose('all');
warning on;
%commment out the following two lines if you want a different trace each
%time you tune the program.
% randn('state',0)
% rand('state',0)

%%--Import and normalize beam shape image
Beamframe = double((imread('BeamShape.tif')));

%imagesc(Beamframe); colormap(gray(256));
figure;
imagesc(Beamframe)
    axis equal tight
    box on
    colormap(gray)
axis xy;
maxVAlue= max(max(Beamframe));

NewBeamImg = Beamframe./maxVAlue;


%%---------This is the file of the data.
fName = 'D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20130529_DnaX_5minFrameRate\Tiffs\FluorescenceRollingBalled.tif';

info = imfinfo(fName);
num_images = numel(info)

%Load first frame
%FrameFirst = imread(fName, 1, 'Info', info);






for k = 1:num_images
    %tStart = tic;
    %fprintf('Image %d is being processed... \n', k);
    %we start here at i=2 since the first frame has been aanlysed already
    fr = imread(fName, k, 'Info', info);
       
    fr = double(imdivide(double(fr),NewBeamImg));
        
    %fNameToWrite = ['/Volumes/LittleMonster/PhD_Docs/20130403_DnaN_MotABDeleted_TimeSeries_37Degrees_S1/TIFFs/FL_RollingBalled_IllumCorrected/FL_RollingBalled_IllumCorrected' num2str(k) '.tif']; 
    fNameToWrite = ['D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20130529_DnaX_5minFrameRate\Tiffs\FL_RollingBalled_IllumCorrected\FL_RollingBalled_IllumCorrected' num2str(k) '.tif']; 
    
    AMIN = 0;
    AMAX = 65535;
    %I = mat2gray(A,[AMIN AMAX])
    
    imwrite(uint16(65535*mat2gray(fr,[AMIN AMAX])),fNameToWrite,'Compression','none');

end







    
    
%     
%     BeamMultiplieddImg = double(immultiply(double(Iout),NewBeamImg));
%     
%     fNameMul = ['Images_Info/NoiseAdded/BeamMultiplied/Img_BeamMultiplied' num2str(i) '.tif']; 
%     imwrite(uint16(65535*mat2gray(BeamMultiplieddImg,[AMIN AMAX])),fNameMul,'Compression','none');




%

% im_sz = [512 512]; %the size of the image
% num_spots =250; % the number of spots in image
% %spot_pos =[64 64]; %positions of where the spots should come
% plotValue = 0;  
% %intensityMu=8.1919e+03;
% intensityMu=2500;
% intensitySTD=intensityMu/10;
% %intensitySTD=0;
% 
% %noiseMuValue = 440; %This value is taken as the one when there is no light on the camera
% noisSTDValue = intensityMu/5;
% noiseMuValue = 0; %This value is taken as the one when there is no light on the camera
% %noisSTDValue = 0;
% 
% 
% widthMu=1.2;
% 
% %%---------These are the  conditions for the trace generator
% samples=200;          % Number of points in the whole trace  
% IntTraces=zeros(num_spots,(samples-1));
% 
% % growthrate=40/25;     % stepsize/window
% 
% Step1=100;         % Size of the steps; 
% noise=0;         % Noise level was 10 Step1/noise is the signal-to-noise
% %SigmaStep1=20;
% window = 15; %was 25, changed it to 15 since in the paper by Sherratt they used a probability of 0.07
% growthrate=Step1/window; % stepsize/window
% nw1=Step1/growthrate;  %what is this? it looks like the "window"
% %nw1=10;
% 
% data=zeros(samples,2);
% data(:,2)=noise*randn(samples,1);  % Use Gaussian noise with amplitude "noise"
% data(:,1)=(1:1:samples)';           % Time axis; this is arbitrary, here time between frames is taken as 255 ms
% %%-------------------------------------------------------------------------
% 
% figure;
% hold on;
% for k=1:num_spots
%     %%-------Make a trace
%     %rng('shuffle')
%     %rng(cputime,'twister')
%     %rand('seed',cputime);
% %stream = RandStream.getGlobalStream
%     data=zeros(samples,2);
% %data(:,2)=noise*randn(samples,1);  % Use Gaussian noise with amplitude "noise"
% %data(:,1)=(1:1:samples)';           % Time axis; this is arbitrary, here time between frames is taken as 255 ms
% 
% %%-------------------------------------------------------------------------
%     stepped = 0;
%         for j=1:samples    
%  
%           if stepped == 0
%             teken=1;  
%             step=Step1*ceil(rand(1,1) - (1-1/nw1));
%     		data(j:samples,2) = data(j:samples,2) + teken*step;
%             if step >0
%                 stepped=1;
%             end;
%           else
%               data(j:samples,2) = data(j:samples,2);
%               %stepped=1;
%           end;
%             %chance of 1/nw that it is 'step'; otherwise 0;
%         end
% 
% signal = data(:,2);
% 
% %%%--- End step maker ---%%%
% 
% %%% Make steps go down
% signal = -signal;
% NewSignal = signal(1:length(signal)-1);
% %%% Add some offset to end around zero
% NewSignal = NewSignal - mean(NewSignal(length(NewSignal)-5:length(NewSignal)));
% NewSignal = NewSignal./max(NewSignal);
% IntTraces(k,:)= NewSignal;
% 
% plot(IntTraces(k,:))
% IntTraces(k,:);
% if (isnan(IntTraces(k,:)))
%     IntTraces(k,:) = 0;
% end;
% 
%     %%-------------------
% end;
% 
% hold off;
% 
% 
% 
% %I = I./(max(I));
% 
% 
% % Xin = [100 150 250 300 350]';
% % Yin = [256 256 256 256 256]';
% % [Iout,Inoise,spot_data]=spotmaker(im_sz,num_spots,'spot_pos',[Xin,Yin],'plot',plotValue,'noise_mu',noiseMuValue,'noise_std',noisSTDValue,'int_mu',intensityMu,'int_std',intensitySTD,'wid_mu',widthMu);
% 
% %[Iout,Inoise,spot_data]=spotmaker(im_sz,num_spots,'plot',plotValue,'noise_mu',noiseMuValue,'noise_std',noisSTDValue,'int_mu',intensityMu,'int_std',intensitySTD,'wid_mu',widthMu);
% 
% %
% 
% for i=1:(samples-1)
%     %intensityMu = intensityMu.*I(i);
%     %Xin = [100 150 250 300 350]';
%     %Yin = [256 256 256 256 256]';
%     %[Iout,Inoise,spot_data]=spotmaker(im_sz,num_spots,'spot_pos',[Xin,Yin],'plot',plotValue,'noise_mu',noiseMuValue,'noise_std',noisSTDValue,'int_mu',intensityMu_array(i),'int_std',intensitySTD,'wid_mu',widthMu);
%    
%     %for k=1:num_spots
%     
%     if ((i==1))
%     [Iout,Inoise,spot_data]=spotmaker_DifferentIntensities(im_sz,num_spots,'plot',plotValue,'noise_mu',noiseMuValue,'noise_std',noisSTDValue,'int_mu',intensityMu,'int_std',intensitySTD,'wid_mu',widthMu,'traceVector', IntTraces(:,i));
%     
%         %[Iout,Inoise,spot_data]=spotmaker(im_sz,num_spots,'plot',plotValue,'noise_mu',noiseMuValue,'noise_std',noisSTDValue,'int_mu',intensityMu,'int_std',intensitySTD,'wid_mu',widthMu,IntTraces(:,i));
%     Xin = spot_data.Xcent';
%     Yin = spot_data.Ycent';
%     else
%        [Iout,Inoise,spot_data]=spotmaker_DifferentIntensities(im_sz,num_spots,'spot_pos',[Xin,Yin],'plot',plotValue,'noise_mu',noiseMuValue,'noise_std',noisSTDValue,'int_mu',intensityMu,'int_std',intensitySTD,'wid_mu',widthMu,'traceVector',IntTraces(:,i));
%        
%      %   [Iout,Inoise,spot_data]=spotmaker(im_sz,num_spots,'spot_pos',[Xin,Yin],'plot',plotValue,'noise_mu',noiseMuValue,'noise_std',noisSTDValue,'int_mu',intensityMu,'int_std',intensitySTD,'wid_mu',widthMu,IntTraces(:,i));
%     end
% 
%     fName = ['Images_Info/NoiseAdded/NoBeamMultilpication/Img_raw' num2str(i) '.tif']; 
%     
%     AMIN = 0;
%     AMAX = 65535;
%     %I = mat2gray(A,[AMIN AMAX])
%     
%     imwrite(uint16(65535*mat2gray(Iout,[AMIN AMAX])),fName,'Compression','none');
%     
%     BeamMultiplieddImg = double(immultiply(double(Iout),NewBeamImg));
%     
%     fNameMul = ['Images_Info/NoiseAdded/BeamMultiplied/Img_BeamMultiplied' num2str(i) '.tif']; 
%     imwrite(uint16(65535*mat2gray(BeamMultiplieddImg,[AMIN AMAX])),fNameMul,'Compression','none');
%     
%     dlmwrite('Images_Info/NoiseAdded/SpotPeakHeightIntensities.txt', [spot_data.ints], 'delimiter', '\t', ...
%         'precision', 6,'-append')
%     dlmwrite('Images_Info/NoiseAdded/SpotPeakWidthIntensities.txt', [spot_data.stds], 'delimiter', '\t', ...
%         'precision', 6,'-append')
% %     dlmwrite('SpothIntensities_InclNoise.txt', [spotIntensitiesNoise], 'delimiter', '\t', ...
% %         'precision', 6,'-append')
%     
%     %imwrite(uint16((Iout)),fName,'Compression','none')
%     %imwrite(uint16(65535*mat2gray(Iout)),fName,'Compression','none')
%     %imwrite(uint16(Iout),fName,'Compression','none')
%     %end;
% end;
% 
% %spot_data.ints(i) = peak height (intensity) of the Gaussian spot i.
% % spot_data.stds(i) = peak width (STD) of the Gaussian spot i.
% 
% 
% dlmwrite('Images_Info/NoiseAdded/SpotPositions.txt', [Xin Yin], 'delimiter', '\t', ...
%         'precision', 6)
% 
% % figure;
% % imagesc(FrameFirst)
% %     axis equal tight
% %     box on
% %     colormap(gray)
% %     axis xy;
%     
%     
%     
%     
% %     figure;
% %     hold on
% %     imagesc(Iout)
% %     axis equal tight
% %     box on
% %     colormap(gray)
% %  
%     
% 


% %Images_IllumCorrector
% clear all
% dipstart
% basepath='D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20130221_HigherFrameRateMeasurement\';
% infiles='Fl_RB\FluRB*';
% outfiles='FL_RBIC\FL_RbIc';
% 
% lb2=strcat(basepath,infiles)
%     ff = readtimeseries(lb2,'tif',[0 405]);              %fluo stack
%    % ff = readtimeseries(lb2,'tif',[1 10]);              
%   
%  lb4=strcat(basepath, 'BeamShape.tif');
%  
%  illum0=readim(lb4,'tif');
%  illum1=illum0/max(illum0(:));
%  
%  [r,c,d]=size(ff);
%  for i=0:d-1
%      i
%      im=squeeze(double(ff(:,:,i)));
%      im_c=uint16(round(im./illum1));
%      nr=num2str(i,'%03.0f');
%      nam=strcat(basepath,outfiles,'Fr',nr,'.tif') ;
%      writeim(im_c,nam,'tif',1);
%      close all
%  end