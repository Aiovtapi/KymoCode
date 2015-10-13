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

%GaussCalcs
%% Histogram Intensities
n=Ncells;
figure(4)
hold on
hist(nonzeros(d{n+1}.IntI),20);
  h = findobj(gca,'Type','patch');
  h.FaceColor = [0 0.7 0.7];
  h.EdgeColor = 'w';
hold off
% axis([0 2100 0 80])

xlabel('Integrated Intensity counts (-)','FontSize',20,'FontWeight','bold');
ylabel('Frequency (-)','FontSize',20,'FontWeight','bold');
title('4th Brightest Spot Frequency of Intensity counts (Tus)','FontSize',20)
txtbox1={TotCellsStr,MeanIntStrD,StdIntStrD};
annotation('textbox', [0.65,0.8,0.1,0.1],...
           'String', txtbox1,'FontSize',20,'FontWeight','bold')
set(gca,'FontSize',20);
   
%%
figure(5)
%hist(cat(1,D{21}.X,nonzeros(D{21}.XR1),nonzeros(D{21}.XR2),nonzeros(D{21}.XR3)),20);
hist(d{n+1}.X,25);
xlabel('X Position (-)','FontSize',20,'FontWeight','bold');
ylabel('Frequency (-)','FontSize',20,'FontWeight','bold');
title('Frequency of spot X Position dif','FontSize',20)
set(gca,'FontSize',20);

%% Plot of position w.r.t. cell
figure(4)
hold on
plot(d{n+1}.X,d{n+1}.Y,'oc');
plot([0 1],[0 0],'r',[0 0],[0,1],'r',[0 1],[1 1],'r',[1 1],[0 1],'r')
hold off
axis([-.1 1.1 -.1 1.1])
xlabel('X Position (-)','FontSize',20);
ylabel('Y Position (-)','FontSize',20);
title('Detected localizations within the cell dif','FontSize',20)

%% Plot w.r.t. life cycle
txtbox7={TotCellsStr};

% figure(7)
% hold on
% plot(xcc,M(:,1),'g',xcc,MR1(:,1),'b',xcc,MR2(:,1),'k',xcc,MR3(:,1),'r');
% hold off
% xlabel('Normalized Cell Time (-)','FontSize',16);
% ylabel('Integrated Intensity (-)','FontSize',16);
% title('Mean Integrated Intensity vs Normalised Cell Time','FontSize',20)
% annotation('textbox', [0.15,0.875,0.1,0.05],...
%            'String', txtbox7,'FontSize',14,'FontWeight','bold')

countsperlabel=990;
figure(8)
hold on
 plot(xcc,(Md(:,1)),'r','LineWidth',4)
plot(xcc,(MdR1(:,1)),':r','LineWidth',2)
plot(xcc,(MdR2(:,1)),':b','LineWidth',2)
% plot(xcc,(MdR2(:,1)),':r','LineWidth',2)
% plot(xccd,MdR3(:,1),'r','LineWidth',4)
% plot(xcc,(M(:,1)),'r','LineWidth',4)
% plot(xccd,Md(:,7)+Md(:,8),'--b')
% plot(xccd,Md(:,7)-Md(:,8),'--b')
% plot(xccd,Md(:,1)+MdR1(:,1)+MdR2(:,1)+MdR3(:,1)+sqrt(Md(:,6).^2+MdR1(:,6).^2+MdR2(:,6).^2+MdR3(:,6).^2),'--r')
% plot(xccd,Md(:,1)+MdR1(:,1)+MdR2(:,1)+MdR3(:,1)-sqrt(Md(:,6).^2+MdR1(:,6).^2+MdR2(:,6).^2+MdR3(:,6).^2),'--r')
hold off
xlabel('Normalized Cell Time (-)','FontSize',20,'FontWeight','bold');
ylabel('Number of Proteins (-)','FontSize',20,'FontWeight','bold');
L=legend('DnaN','std');
set(L,'FontSize',20);
set(gca,'FontSize',20);
title('Number of Bound Proteins vs Normalised Cell Time (DnaN)','FontSize',24)
annotation('textbox', [0.135,0.88,0.08,0.03],...
           'String', txtbox7,'FontSize',20,'FontWeight','bold')
       
%% Tus w.r.t. dif

Fac=1/(length(xccd)/(MeanCellLifed));

figure(8)
hold on
% plot(M(:,2),xcc/Fac,'b','LineWidth',4)%M(:,2)-M(:,3),xcc,'b',M(:,2)+M(:,3),xcc,'b',...
%    plot(MR1(:,2),xcc/Fac,'b','LineWidth',4)
   plot(MR2(:,2),xcc/Fac,'b','LineWidth',4)
%   plot(MR3(:,2),xcc/Fac,'b','LineWidth',4)
 plot(Md(:,2),xccd/Fac,'r','LineWidth',4)%,Md(:,2)-Md(:,3),xccd,'r',Md(:,2)+Md(:,3),xccd,'r')
%  plot(MdR1(:,2),xccd/Fac,'r','LineWidth',4)
%   plot(MdR2(:,2),xccd/Fac,'r','LineWidth',4)
%     plot(MdR3(:,2),xccd/Fac,'r','LineWidth',4)
line([0 1],[1 1],'LineWidth',2)
line([0.5 0.5],[1 2])
% hline2=refline([0 MeanCellLifeT*5]);
% hline2.Color='k';
axis([0 1 0 2])
hold off

ylabel('Normalized Cell Time (-)','FontSize',20);
xlabel('X Position (-)','FontSize',20);
title('Mean Normalised X position vs Normalised Cell Time','FontSize',20)
legend('Tus','dif','Div time')
set(legend,'FontSize',18)

%% Ratio of activity
figure(9)
hold on
plot(xcc,Mdratio,'LineWidth',4);
% plot(xccd,Mdratio,'LineWidth',4);
hold off
xlabel('Normalized Cell Time (-)','FontSize',20,'FontWeight','bold');
ylabel('Ratio (-)','FontSize',20,'FontWeight','bold');
L=legend('Tus','dif');
set(L,'FontSize',20);
set(gca,'FontSize',20);
title('Time dependence - ratio of spot vs cell fluorescence','FontSize',24)
annotation('textbox', [0.6,0.85,0.1,0.05],...
           'String', txtbox7,'FontSize',20,'FontWeight','bold')

%% Another
figure(10)
for i=2
hold on
plot((1:length(d{i}.x(:,6)))/length(d{i}.x(:,6)),d{i}.x(:,6)+d{i}.xR1(:,6)+d{i}.xR2(:,6),'b')
% plot((1:length(D{i}.x(:,6)))/length(D{i}.x(:,6)),D{i}.x(:,6),'r');
hold off
end
