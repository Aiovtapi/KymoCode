close all
clear Sim RadDiff SpotWriter

%% load data
TestFolder='GaussNoiseTest4';
ImageFolder='images220nmPS/';

for N=13:29
    
clear DistanceNM MatIndex Distance SNR_Sens_PPV RadX RadSim RadDiff
LionFitFile=num2str(N);
SaveFile=LionFitFile;

NoiseData=load(strcat('/Users/rleeuw/Work/DataAnalysis/BlurLab/',TestFolder,'/',TestFolder,'-2.txt'));

load(strcat('/Users/rleeuw/Work/DataAnalysis/BlurLab/',TestFolder,'/',ImageFolder,'/Results/',LionFitFile));

Sim=NoiseData;

% Blur parameters 

Border=500;
nmperpixel=220;

OffsetX=149;
OffsetY=501.5;

% Blur Box Dimensions (?m)
[YSize,XSize]=size(ydatacrpd{1,1});

Lx=XSize*nmperpixel/1000;
Ly=YSize*nmperpixel/1000;
Lz=0.05;

ShiftX=OffsetX/nmperpixel;
ShiftY=OffsetY/nmperpixel;

frame = 1;

%Normalize:
Sim(:,1)=Sim(:,1)/Lx;
Sim(:,2)=Sim(:,2)/Ly;

%Shift accordingly:
Sim(:,1)=Sim(:,1)*XSize+ShiftX;
Sim(:,2)=Sim(:,2)*YSize+ShiftY;

%calculate distances between spots and simulation points

for j=1:NSpots
RadX(j,1)=sqrt(x{j}(1,2)^2+x{j}(1,4)^2);
end
RadSim(:,1)=sqrt(Sim(:,1).^2+Sim(:,2).^2);

for i=1:NSpots
    for j=1:length(Sim(:,1))   
        RadDiff(i,j)=abs(RadX(i,1)-RadSim(j,1));
    end    

end
    [MatIndex,Distance]=LionLAP(RadDiff);

    DistanceNM=Distance*80;

%The PSF is 80x80 nm. I will use criterium that if detection is more 120 nm
%away, it is a FALSE positive.

DistanceNM(DistanceNM>160)=NaN;
%DistanceNM now contains True positives.
%Number of True negatives 

Sensitivity=(sum(~isnan(DistanceNM)))/(length(Sim(:,1)));
PPV=sum(~isnan(DistanceNM))/(NSpots);

SNR_Sens_PPV=[mean(SNR) Sensitivity PPV];

% figure(1)
% hold on
% plot(sqrt(Sim(:,1).^2+Sim(:,2).^2),'r','LineWidth',3);
% plot(sqrt(x{1}(:,2).^2+x{1}(:,4).^2),'b','LineWidth',3);
% hold off
% legend('Simulation','LionFit');
% title('Single spot with S/N = 1000','FontSize',20);
% xlabel('Frame (-)','FontSize',16);
% ylabel('Radial position from origin (-)','FontSize',16);
% 
% 
% figure(2)
% hold on
% imagesc(ydatacrpd{frame,1});
% plot(Sim(frame,1),Sim(frame,2),'or','LineWidth',3);
% plot(x{1}(frame,2),x{1}(frame,4),'+k','LineWidth',3);
% hold off
% legend('Simulation','LionFit','FontSize',16);
% title('Comparison of localizations of sim and fit','FontSize',16);
% 
% figure(3)
% plot(SNRBlur,'LineWidth',3);
% xlabel('Frame (-)','FontSize',16)
% ylabel('SNR (-)','FontSize',16)
% legend('SNR','FontSize',16)
% title('SNR theoretical during the simulation','FontSize',20)
% 
figure(3)
hold on
imagesc(ydatacrpd{1});
plot(Sim(:,1),Sim(:,2),'+r','MarkerSize',10)
for j=1:NSpots
plot(x{j}(1,2),x{j}(1,4),'+w','MarkerSize',10)
end
hold off
save(strcat('SimulationResults/KAKTEST/',SaveFile),'SNR_Sens_PPV');
end
