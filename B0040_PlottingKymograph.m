clc;
clear all;

warning on;


exp='001_DnaN_TUS_dif_30122014_DnaNsignal';
initval=A001_Images_Set_Experiment(exp);

%KymoNr = 1;

SelectedKymoGraphName = strcat('Kymograph_FLExp001_DnaN_TUS_dif_30122014_DnaNsignalChan_x245')%,num2str(KymoNr));


Kymo=double((imread(strcat(initval.basepath,'Kymographs/',SelectedKymoGraphName,'.tif')))); %kymograph


%%--Import and normalize beam profile image
%Beamframe = double((imread('BeamShape.tif')));

H = figure;
%imagesc(Kymo)
%Kymo = flipud(Kymo);
pcolor(Kymo)
axis equal tight
%box on
colormap(hot)
shading flat;
axis xy;
set(gca,'YDir','reverse');


   set(gca, 'fontsize', 26, 'linewidth', 2, 'fontweight', 'bold','box','off');
   title('Fluorescence Kymograph dif site', 'fontsize', 18, 'fontweight', 'bold')
   ylabel('Time (min)', 'fontsize', 26, 'fontweight', 'bold')
   xlabel('Position (px)', 'fontsize', 26, 'fontweight', 'bold')
   %a = annotation('textbox','Position',[0.65 0.65 0.35 0.16], 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('N = %d',N_CombinedSeries) );
   %axis([0 100 0 55e-3])
   %set(gca, 'XTick', [-2e3 0 2e3 4e3 6e3 8e3 10e3], 'YTick', [0 10e-5 20e-5 30e-5]); 
   set(gca, 'XTick', [0 20 60 100 140]);
   %set(gca,'TickLength',[0.02 0.02]);
   SavingPathAndName=strcat(initval.basepath,'Kymographs\','CFP_Fluorescence_Kymograph1')% ,num2str(KymoNr)); %kymograph
   print(H, '-dpng','-zbuffer', '-r600',SavingPathAndName)

%kymo_BF=fliplr(kymo_BF);
%kymo_FL=fliplr(kymo_FL);
%figure; pcolor(kymo_BF); shading flat; colormap hot;
%figure; pcolor(kymo_FL); shading flat; colormap hot;


% NewBeamImg = Beamframe./maxVAlue;
% 
% %This is the fluorescence stack to be corrected
% fName ='/data/FluorescenceStackToBeCorrected.tif';
% 
% info = imfinfo(fName);
% num_images = numel(info)
% 
% 
% for k = 1:num_images
%     fr = imread(fName, k, 'Info', info);
%        
%     fr = double(imdivide(double(fr),NewBeamImg));
%         
%     fNameToWrite = ['/IllumCorrectedFluorescenceImages/FL_Corrected_Images' num2str(k) '.tif']; 
%     
%     AMIN = 0;
%     AMAX = 65535;
%     
%     imwrite(uint16(65535*mat2gray(fr,[AMIN AMAX])),fNameToWrite,'Compression','none');
% 
% end