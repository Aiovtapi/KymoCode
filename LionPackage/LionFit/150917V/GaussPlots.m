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
% IntensitySumSd=0;

i=1;

Resd=cell(length(dXmat{1}(:,1)),1);
ResT=cell(length(TXmat{1}(:,1)),1);
TotalLengthResd=0;
TotalLengthResT=0;

for k=1:length(dXmat{1}(:,1))
Resd{k}=dXmat{i}(k,IdXmat{k,i});
TotalLengthResd=TotalLengthResd+size(Resd{k},2);
end

for k=1:length(TXmat{1}(:,1))
ResT{k}=TXmat{i}(k,ITXmat{k,i});
TotalLengthResT=TotalLengthResT+size(ResT{k},2);
end


Time_Resd=(1:size(Resd,1));
StretchedResd=cell(max(SpotNumbersd),1);
delta=0;
Time_Resd=Time_Resd';
skip=0;

Time_ResT=(1:size(ResT,1));
StretchedResT=cell(max(SpotNumbersT),1);
deltaT=0;
Time_ResT=Time_ResT';
skipT=0;

% for k=1:size(Resd,1)
%    
%     L=size(Resd{k},2);
%     
%     if L>1
%         for p=1:L   
%         StretchedResd{k+delta+p-1)=Resd{k}(p);
%         end
%         for n=1:p-1
%         Time_Resd=[Time_Resd(1:k+delta);Time_Resd(k+delta);Time_Resd(k+1+delta:end)];
%         end
%         delta=delta+L-1; %by chance this is the occurances of multiple spots. This is how many times we shift the index for an occurance of multiple spots.
%     else
%         StretchedResd(k+delta)=Resd{k};
%     end
% end

for k=1:size(Resd,1)
   
    L=size(Resd{k},2);
    
    if L>1
        for p=1:L
        StretchedResd{p}(k,2)=Resd{k}(p);
        StretchedResd{p}(k,1)=Time_Resd(k);
        end
    else
        StretchedResd{1}(k,1)=Time_Resd(k);
        StretchedResd{1}(k,2)=Resd{k}(1);
    end
end

for k=1:size(ResT,1)
   
    L=size(ResT{k},2);
    
    if L>1
        for p=1:L
        StretchedResT{p}(k,2)=ResT{k}(p);
        StretchedResT{p}(k,1)=Time_ResT(k);
        end
    else
        StretchedResT{1}(k,1)=Time_ResT(k);
        StretchedResT{1}(k,2)=ResT{k}(1);
    end
end


 for n=1:max(SpotNumbersd)

     [i,j,v]=find(StretchedResd{n});
      
     LNZ=length(nonzeros(StretchedResd{n}))/2;
     
     StretchedResd{n}=[];
     
     I=linspace(1,LNZ,LNZ)';
     I=[I;I];

     J=ones(LNZ,1);
     J=[J;ones(LNZ,1)*2];
     
     for m=1:length(i)
     StretchedResd{n}(I(m),J(m))=v(m);
     end

     
 end

 for n=1:max(SpotNumbersT)

     [i,j,v]=find(StretchedResT{n});
      
     LNZ=length(nonzeros(StretchedResT{n}))/2;
     
     StretchedResT{n}=[];
     
     I=linspace(1,LNZ,LNZ)';
     I=[I;I];

     J=ones(LNZ,1);
     J=[J;ones(LNZ,1)*2];
     
     for m=1:length(i)
     StretchedResT{n}(I(m),J(m))=v(m);
     end

     
 end

figure(1)
hold on
for K=1
scatter(StretchedResd{K}(:,1)/(length(StretchedResd{K}(:,1))),StretchedResd{K}(:,2),[],'filled')
end
for K=1:3
scatter(StretchedResT{K}(:,1)/(length(StretchedResT{K}(:,1))),StretchedResT{K}(:,2),[],'filled')
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
scatter(0:1/(MeanBacLifeT-4):1,T{i}.XNorm{j}(:,2),T{i}.x{j}(:,1)/10,'or','filled');
scatter(0:1/(MeanBacLifeT-4):1,S{i}.x{j}(:,2),S{i}.x{j}(:,1)/40,'ob','filled');
end
hold off
axis([0 1 0 1])
xlabel('Normalized Cell Time (-)');
ylabel('Normalized Cell Position (-)');
set(gca,'FontSize',20);

%% Distance between tus and dif vs. time

i=1;

for n=1:Ncells
    DistanceMeanMatrixAvgpFrame=nanmean(ddweighted{n},2);
    DistanceStdMatrix=nanstd(ddweighted{n},1,2);
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
plot(0:2/(MeanBacLifeT+3):2,S{3}.x{3}(:,2),'b','LineWidth',5);
plot(0:2/(MeanBacLifeT+3):2,Sd{3}.x{1}(:,2),'r','LineWidth',5);
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

