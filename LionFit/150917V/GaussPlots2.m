%% x position of spots
close all

dif.Nspots=size(dif.I,2);

figure(1)

hold on
scatter(0:1/(MeanBacLifed-1):1,dif.x{1},dif.I{1}/50+1,'b','filled')
% scatter(0:1/(MeanBacLifed-1):1,dif.x{2},dif.I{2}/50+1,'r','filled')
% scatter(0:1/(MeanBacLifed-1):1,dif.x{3},dif.I{3}/50+1,'k','filled')
% plot(0:1/(MeanBacLifed-1):1,dif.x{1}+dif.xstd{1},'b','LineWidth',1)
% % plot(0:1/(MeanBacLifed-1):1,dif.x{1}-dif.xstd{1},'b','LineWidth',1)
scatter(0:1/(MeanBacLifeT-1):1,Tus.x{1},Tus.I{1}/20+1,'r','filled')
scatter(0:1/(MeanBacLifeT-1):1,Tus.x{2},Tus.I{2}/20+1,'m','filled')
scatter(0:1/(MeanBacLifeT-1):1,Tus.x{3},Tus.I{3}/20+1,'k','filled')
% scatter(0:1/(MeanBacLifeT-1):1,Tus.x{4},Tus.I{4}/20+1,'g','filled')
% plot(0:1/(MeanBacLifed-1):1,Tus.xmean+Tus.xmeanstd,'r','LineWidth',1)
% plot(0:1/(MeanBacLifed-1):1,Tus.xmean-Tus.xmeanstd,'r','LineWidth',1)
% scatter(0:1/(MeanBacLifed-1):1,dif.x{2},dif.I{2}/50+1,'b','LineWidth',2)
% scatter(0:1/(MeanBacLifed-1):1,dif.x{3},dif.I{3}/50+1,'b','LineWidth',2)
% scatter(0:1/(MeanBacLifeT-1):1,Tus.x{1},Tus.I{1}/100+1,'r','LineWidth',2)
% scatter(0:1/(MeanBacLifeT-1):1,Tus.x{2},Tus.I{2}/100+1,'r','LineWidth',2)
% scatter(0:1/(MeanBacLifeT-1):1,Tus.x{3},Tus.I{3}/100+1,'r','LineWidth',2)
% scatter(0:1/(MeanBacLifeT-1):1,Tus.x{4},Tus.I{4}/100+1,'r','LineWidth',2)
%  scatter(0:1/(MeanBacLifeT-1):1,Tus.x{5},Tus.I{5}/100+1,'r','LineWidth',2)
%  scatter(0:1/(MeanBacLifeT-1):1,Tus.x{6},Tus.I{6}/100+1,'r','LineWidth',2)
%   scatter(0:1/(MeanBacLifeT-1):1,Tus.x{5},Tus.I{7}/100+1,'r','LineWidth',2)
%  scatter(0:1/(MeanBacLifeT-1):1,Tus.x{6},Tus.I{8}/100+1,'r','LineWidth',2)
hold off
axis([0 1 0.01 1])
xlabel('Normalised cell time (-)','FontSize',16)
ylabel('Normalised long axis position (-)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
legend('dif','Tus','Tus','Tus')
title('Cellular position of Tus proteins and the dif loci')

%% Intensities

%calibration value Tus (by Sriram - see Tus chapter)
Pv=978.7;

figure(2)
hold on
plot(0:1/(MeanBacLifeT-1):1,Tus.FC{1}/Pv,'r','LineWidth',5)
plot(0:1/(MeanBacLifeT-1):1,(Tus.FC{1}+Tus.FCstd{1})/Pv,'r','LineWidth',1)
plot(0:1/(MeanBacLifeT-1):1,(Tus.FC{1}-Tus.FCstd{1})/Pv,'r','LineWidth',1)
% scatter(0:1/(MeanBacLifed-1):1,Tus.I{2},'r','LineWidth',2)
hold off
xlabel('Normalised cell time (-)','FontSize',16)
ylabel('Tus proteins in cell (-)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
title('Tus proteins in the cell vs. time')

%% Tus proteins in a spot:
figure(3)
hold on
plot(0:1/(MeanBacLifeT-1):1,(Tus.I{1}/Pv),'r','LineWidth',5)
plot(0:1/(MeanBacLifeT-1):1,(Tus.I{1}-Tus.Istd{1})/Pv,'r','LineWidth',1)
plot(0:1/(MeanBacLifeT-1):1,(Tus.I{1}+Tus.Istd{1})/Pv,'r','LineWidth',1)
hold off
axis([0 1 0 3])
xlabel('Normalised cell time (-)','FontSize',16)
ylabel('Tus proteins in spot (-)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
title('Tus proteins in a spot vs. time')

%% Total active Tus spots:
figure(4)
hold on
scatter(0:1/(MeanBacLifeT-1):1,(Tus.activespots),100,'b','filled')
hold off
axis([0 1 0.01 10])
xlabel('Normalised cell time (-)','FontSize',16)
ylabel('Number of spots (-)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
title('Number of spots vs. time')

%% Active fraction of Tus

Tus.IoverFCstd=(1/sqrt(3).*(Tus.I{1}./Tus.FC{1})).*sqrt((Tus.Istd{1}./(Tus.I{1})).^2+ ...
    (Tus.FCstd{1}./(Tus.FC{1})).^2);

figure(5)
hold on
plot(0:1/(MeanBacLifeT-1):1,(Tus.I{1}./Tus.FC{1})*100,'b','LineWidth',3)
plot(0:1/(MeanBacLifeT-1):1,((Tus.I{1}./Tus.FC{1})+Tus.IoverFCstd)*100,'b','LineWidth',1)
plot(0:1/(MeanBacLifeT-1):1,((Tus.I{1}./Tus.FC{1})-Tus.IoverFCstd)*100,'b','LineWidth',1)
hold off
axis([0 1 0.01 40])
xlabel('Normalised cell time (-)','FontSize',16)
ylabel('Active fraction (%)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
title('Active fraction in single spot vs. time')

%% Active fraction of Tus

Tus.IoverFCstd=(1/sqrt(3).*(Tus.I{1}./Tus.FC{1})).*sqrt((Tus.Istd{1}./(Tus.I{1})).^2+ ...
    (Tus.FCstd{1}./(Tus.FC{1})).^2);

figure(6)
hold on
plot(0:1/(MeanBacLifeT-1):1,(Tus.I{1}+Tus.I{2}+Tus.I{3}+Tus.I{4}./Tus.FC{1}),'b','LineWidth',3)
plot(0:1/(MeanBacLifeT-1):1,((Tus.I{1}+Tus.I{2}+Tus.I{3}+Tus.I{4}./Tus.FC{1})+Tus.IoverFCstd)*100,'b','LineWidth',1)
plot(0:1/(MeanBacLifeT-1):1,((Tus.I{1}+Tus.I{2}+Tus.I{3}+Tus.I{4}./Tus.FC{1})-Tus.IoverFCstd)*100,'b','LineWidth',1)
hold off
axis([0 1 0.01 6000])
xlabel('Normalised cell time (-)','FontSize',16)
ylabel('Active fraction (%)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
title('Active fraction in single spot vs. time')

%% individual growth curves 
figure(7)
hold on
for i=1:Ncells
plot((1:length(dif.length{i}))*2.5,dif.length{i}/(dif.length{i}(end)),'LineWidth',2)
end
axis([0 120 0.3 1])
xlabel('Time (min)','FontSize',16)
ylabel('Relative growth (-)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
title('Growth curves of relative growth vs. time')
%% genertion hist
for i=1:Ncells
    GenerationTime(i)=size(Sd{i}.x{1},1)*2.5;
end
hist(GenerationTime,10);
axis([0 500 0.01 10])
xlabel('Time (min)','FontSize',16)
ylabel('Frequency (-)','FontSize',16)
set(gca,'FontSize',16,'FontWeight','bold')
title('Generation Time Distribution')

%% Test

hold on
plot(0:1/(length(d{1}.x{1}(:,2))-1):1,d{1}.x{1}(:,2),'b');
plot(0:1/(length(d{1}.x{2}(:,2))-1):1,d{1}.x{2}(:,2),'b');
plot(0:1/(length(d{1}.x{3}(:,2))-1):1,d{1}.x{3}(:,2),'b');
plot(0:1/(length(T{1}.x{1}(:,2))-1):1,T{1}.x{1}(:,2),'r');
plot(0:1/(length(T{1}.x{2}(:,2))-1):1,T{1}.x{2}(:,2),'r');
plot(0:1/(length(T{1}.x{3}(:,2))-1):1,T{1}.x{3}(:,2),'r');