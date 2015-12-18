%%              Overview of Variables
% ------------------------------------------------------------------------
% D{1-n} = cellular data [Amp xpos sigx ypos sigy] (pos normalized to cell)
% "R1,R2,R3" addition to var = 2nd,3rd,4th spot respectively.
% "d" addition to var = dif spot results.
% D{n+1}.X = x position first (brightest) spot. 
% D{n+1}.Y = y position first (brightest) spot.
% D{n+1}.IntI = Integrated Intensity.
% M = Mean value at time points 
% M = [IntSpot Xpos stdX Ypos stdY stdIntSpot IntCell stdIntCell]

% GaussCalcs2

%% dif Position vs. time

i=15;

figure(1)
hold on
for j=1:10
scatter(0:1/(MeanBacLifed+7):1,d{i}.XNorm{j}(:,2),d{i}.x{j}(:,8)/100,'or','LineWidth',3);
scatter(0:1/(MeanBacLifed-1):1,Sd2{i}.x{j}(:,2),Sd2{i}.x{j}(:,1)/100,'ob','LineWidth',3);
end
hold off
axis([0 1 0 1])
xlabel('Normalized Cell Time (-)');
ylabel('Normalized Cell Position (-)');
set(gca,'FontSize',20);

%% Tus Position vs. time

i=15;

figure(1)
hold on
for j=1:10
scatter(0:1/(MeanBacLifeT+2):1,T{i}.XNorm{j}(:,2),T{i}.x{j}(:,1)/10,'or','LineWidth',3);
scatter(0:1/(MeanBacLifeT-1):1,S{i}.x{j}(:,2),S{i}.x{j}(:,1)/40,'ob','LineWidth',3);
end
hold off
axis([0 1 0 1])
xlabel('Normalized Cell Time (-)');
ylabel('Normalized Cell Position (-)');
set(gca,'FontSize',20);

%% Distance between tus and dif vs. time

i=1;

for n=1:Ncells
    DistanceMeanMatrixAvgpFrame=nanmean(ddweighted,2);
    DistanceStdMatrix=nanstd(ddweighted,1,2);
end

figure(3)
hold on
for j=1:10
plot(0:1/(length(DistanceMeanMatrixAvgpFrame)-1):1,DistanceMeanMatrixAvgpFrame,'-k','LineWidth',2);
plot(0:1/(length(DistanceMeanMatrixAvgpFrame)-1):1,DistanceMeanMatrixAvgpFrame+DistanceStdMatrix,'--b','LineWidth',2);
plot(0:1/(length(DistanceMeanMatrixAvgpFrame)-1):1,DistanceMeanMatrixAvgpFrame-DistanceStdMatrix,'--b','LineWidth',2);
end
hold off
% axis([0 1 0 1])
xlabel('Normalized Cell Time (-)');
ylabel('Normalized Cell Position (-)');
set(gca,'FontSize',20);

%% Tus VS dif

figure(3)
hold on
plot(0:2/(MeanBacLifeT-1):2,S{3}.x{3}(:,2),'b','LineWidth',5);
plot(0:2/(MeanBacLifed-1):2,Sd{3}.x{1}(:,2),'r','LineWidth',5);
hold off
axis([0 2 0 1])
xlabel('Normalized Cell Time (-)');
ylabel('Normalized Cell Position (-)');
set(gca,'FontSize',20);
legend('Tus','dif');

%% Histogram Intensities Tus
n=Ncells;

HisIntT=[];

for i=1:n
for j=1:Nspots
    HisIntT=[HisIntT S{i}.x{j}(:,1)];
end
end

figure(4)
hold on
hist(nonzeros(HisIntT),20);
  h = findobj(gca,'Type','patch');
  h.FaceColor = [1 0 0];
  h.EdgeColor = 'w';
hold off

xlabel('Integrated Intensity counts (-)','FontSize',20,'FontWeight','bold');
ylabel('Frequency (-)','FontSize',20,'FontWeight','bold');
title('Spot Intensity counts (Tus)','FontSize',20)
txtbox1={TotCellsStr,sprintf('Mean = %g',nanmean(HisIntT(:))), ...
    sprintf('Std = %g',nanstd(HisIntT(:))),sprintf('Npoints = %g', NpointsT)};
annotation('textbox', [0.65,0.8,0.1,0.1],...
           'String', txtbox1,'FontSize',20,'FontWeight','bold')
set(gca,'FontSize',20);

%% Histogram Intensities dif
n=Ncells;

HisIntd=[];
for i=1:n
for j=1:Nspots
    HisIntd=[HisIntd Sd2{i}.x{j}(:,1)];
end
end

figure(3)
hold on
hist(nonzeros(HisIntd),20);
  h = findobj(gca,'Type','patch');
  h.FaceColor = [0 0 1];
  h.EdgeColor = 'w';
hold off

xlabel('Integrated Intensity counts (-)','FontSize',20,'FontWeight','bold');
ylabel('Frequency (-)','FontSize',20,'FontWeight','bold');
title('Spot Intensity counts (dif)','FontSize',20)
txtbox1={TotCellsStr,sprintf('Mean = %g',nanmean(HisIntd(:))), ...
    sprintf('Std = %g',nanstd(HisIntd(:))),sprintf('Npoints = %g', NpointsT)};
annotation('textbox', [0.65,0.8,0.1,0.1],...
           'String', txtbox1,'FontSize',20,'FontWeight','bold')
set(gca,'FontSize',20);

%%
figure(5)

Weight=cell(Nspots,1);

HisXd=[];

for i=1:Ncells
    for j=1:Nspots
        
        WeightedAverage(i,j)=sum(Sd{i}.x{j}(:,2).*Sd{i}.x{j}(:,1))/sum(Sd{i}.x{j}(:,1));
        
        HisXd=[HisXd d{n+1}.X{j}];
    
    end
end


hist(WeightedAverage(:),9);
xlabel('X Position (-)','FontSize',20,'FontWeight','bold');
ylabel('Frequency (-)','FontSize',20,'FontWeight','bold');
title('Weighted Spot X Position (dif)','FontSize',20)
axis([0 1 0 21]);
set(gca,'FontSize',20);

%% Plot of position w.r.t. cell
j=10;
figure(4)
hold on
plot(d{n+1}.X{j},d{n+1}.Y{j},'oc');
plot([0 1],[0 0],'r',[0 0],[0,1],'r',[0 1],[1 1],'r',[1 1],[0 1],'r')
hold off
axis([-.1 1.1 -.1 1.1])
xlabel('X Position (-)','FontSize',20);
ylabel('Y Position (-)','FontSize',20);
title('Detected localizations within the cell dif','FontSize',20)

%% Plot w.r.t. life cycle
txtbox7={TotCellsStr};
xcc=linspace(0,1,MeanBacLifeT);
xccd=linspace(0,1,MeanBacLifed+1);

% figure(7)
% hold on
% plot(xcc,M(:,1),'g',xcc,MR1(:,1),'b',xcc,MR2(:,1),'k',xcc,MR3(:,1),'r');
% hold off
% xlabel('Normalized Cell Time (-)','FontSize',16);
% ylabel('Integrated Intensity (-)','FontSize',16);
% title('Mean Integrated Intensity vs Normalised Cell Time','FontSize',20)
% annotation('textbox', [0.15,0.875,0.1,0.05],...
%            'String', txtbox7,'FontSize',14,'FontWeight','bold')
j=1;
figure(9)
hold on
plot(0:(1/(length(M{j}(:,1))-1)):1,(M{j}(:,1))/1100,'b','LineWidth',4)
plot(0:(1/(length(M{j}(:,1))-1)):1,M{j}(:,7)/1100,'r','LineWidth',4)
% plot(xcc,M{j}(:,7)-M{j}(:,8),'--b')
% scatter(xcc,(M{j}(:,1)+M{j+1}(:,1)+M{j+2}(:,1)+M{j+3}(:,1)+M{j+4}(:,1)+M{j+5}(:,1)...
%     +M{j+6}(:,1)+M{j+7}(:,1)+M{j+8}(:,1)+M{j+9}(:,1))/1100,'r','LineWidth',4)
% plot(xcc,M{j}(:,1)/200,'r','LineWidth',4)
% plot(xcc,(M{j}(:,1)+M{j+1}(:,1)+M{j+2}(:,1))./M{j}(:,7),'r','LineWidth',4)
hold off
xlabel('Normalized Cell Time (-)','FontSize',20,'FontWeight','bold');
ylabel('Number of Tus Protein (-)','FontSize',20,'FontWeight','bold');
L=legend('Spots','Total');
set(L,'FontSize',20);
set(gca,'FontSize',20);
title('Number of Tus Protein vs Normalised Cell Time (dif)','FontSize',24)
annotation('textbox', [0.25,0.85,0.1,0.05],...
           'String', txtbox7,'FontSize',20,'FontWeight','bold')
       
%% Tus w.r.t. dif
txtbox7={TotCellsStr};
xcc=linspace(0,1,MeanBacLifeT);
xccd=linspace(0,1,MeanBacLifed);

Fac=1/(length(xccd)/(MeanCellLifed+10));

SPOTTESTER=[];
for j=1:Nspots
SPOTTESTER=[SPOTTESTER Md{j}(:,1)];
end


figure(8)
hold on
% plot(M(:,2),xcc/Fac,'b','LineWidth',4)%M(:,2)-M(:,3),xcc,'b',M(:,2)+M(:,3),xcc,'b',...
%    plot(MR1(:,2),xcc/Fac,'b','LineWidth',4)
%   plot(MR2(:,2),xcc/Fac,'b','LineWidth',4)
%   plot(MR3(:,2),xcc/Fac,'b','LineWidth',4)

 plot(xccd,nanmean(SPOTTESTER,2),'b','LineWidth',4)

%  plot(MdR1(:,2),xccd/Fac,'r','LineWidth',4)
%   plot(MdR2(:,2),xccd/Fac,'r','LineWidth',4)
%     plot(MdR3(:,2),xccd/Fac,'r','LineWidth',4)
% line([0 1],[1 1],'LineWidth',2)
% line([0.5 0.5],[1 2])
% hline2=refline([0 MeanCellLifeT*5]);
% hline2.Color='k';
% axis([0 1 0 2])
hold off

ylabel('Integrated Intensity (-)','FontSize',20);
xlabel('Normalized Cell Time (-)','FontSize',20);
title('Mean Spot Intensity vs Normalised Cell Time','FontSize',20)
legend('dif','dif','Div time')
set(legend,'FontSize',18)

%% Ratio of activity
figure(9)
hold on
% plot(xcc,Mratio,'LineWidth',4);
plot(xccd,Md(:,1),'LineWidth',4);
hold off
xlabel('Normalized Cell Time (-)','FontSize',20,'FontWeight','bold');
ylabel('Ratio (-)','FontSize',20,'FontWeight','bold');
L=legend('Tus','dif');
set(L,'FontSize',20);
set(gca,'FontSize',20);
title('Time dependence - ratio of spot vs cell fluorescence','FontSize',24)
annotation('textbox', [0.6,0.85,0.1,0.05],...
           'String', txtbox7,'FontSize',20,'FontWeight','bold')

