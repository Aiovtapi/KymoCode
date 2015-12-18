%% DnaN filtering 

GaussCalcs
%% 

LB=8000;
countsperlabel=990;

FKi=zeros(MaxBacLife,Ncells);
FKiR1=zeros(MaxBacLife,Ncells);
FKiR2=zeros(MaxBacLife,Ncells);


for i=1:Ncells
    DumKi=Ki(:,i)-LB;
    DumKi(DumKi<0)=0;
    DumKi(DumKi>0)=DumKi(DumKi>0)+LB;
    FKi(:,i)=DumKi;
    
    DumKiR1=KiR1(:,i)-LB;
    DumKiR1(DumKiR1<0)=0;
    DumKiR1(DumKiR1>0)=DumKiR1(DumKiR1>0)+LB;
    FKiR1(:,i)=DumKiR1;
    
    DumKiR2=KiR2(:,i)-LB;
    DumKiR2(DumKiR2<0)=0;
    DumKiR2(DumKiR2>0)=DumKiR2(DumKiR2>0)+LB;
    FKiR2(:,i)=DumKiR2;
    
end

for j=1:MaxBacLife
    FM(j,1)=mean(FKi(j,:));
    FMR1(j,1)=mean(FKiR1(j,:));
    FMR2(j,1)=mean(FKiR2(j,:));
    
    FM(j,6)=std(FKi(j,:));
    FMR1(j,6)=std(FKiR1(j,:));
    FMR2(j,6)=std(FKiR2(j,:));
end

%% Individual Intensity plots
% plots 2 8 10 12 13

for i=2
figure(i)
hold on
plot(xcc,FKi(:,i),'ob','LineWidth',4);
plot(xcc,FKiR1(:,i),'or','LineWidth',4);
plot(xcc,FKi(:,i)+FKiR1(:,i),'ok','LineWidth',4);
hold off
end

%% Overall Intensity Plots
txtbox7={TotCellsStr};
figure(2)
hold on
plot(xcc,(FM(:,1)+FMR1(:,1))/countsperlabel,'r','LineWidth',4)
plot(xcc,(FM(:,1)+FMR1(:,1)+FM(:,6))/countsperlabel,':r','LineWidth',1)
plot(xcc,(FM(:,1)+FMR1(:,1)-FM(:,6))/countsperlabel,':r','LineWidth',1)
hold off
xlabel('Normalized Cell Time (-)','FontSize',20,'FontWeight','bold');
ylabel('Number of Proteins (-)','FontSize',20,'FontWeight','bold');
L=legend('DnaN','std');
set(L,'FontSize',20);
set(gca,'FontSize',20);
title('Number of Bound Proteins vs Normalised Cell Time (DnaN)','FontSize',24)
annotation('textbox', [0.135,0.88,0.08,0.03],...
           'String', txtbox7,'FontSize',20,'FontWeight','bold');
%% DNAN
for i=1
figure(i)
hold on
% plot(xcc,Ki(:,i),'ob','LineWidth',4);
% plot(xcc,KiR1(:,i),'or','LineWidth',4);
plot(xcc,FM(:,i)+FM(:,i),'ok','LineWidth',4);
hold off
end