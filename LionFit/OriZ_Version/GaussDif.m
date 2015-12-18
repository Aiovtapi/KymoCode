%% Setting criteria for dif spots

GaussCalcs
%% 

LB=5000;

FKdi=zeros(MaxBacLifed,Ncells);
FKdiR1=zeros(MaxBacLifed,Ncells);
FKdiR2=zeros(MaxBacLifed,Ncells);
MK=zeros(MaxBacLifed,Ncells);

Dcrit=0.2; % minimal distance between spots for identification

% Filtering -- lower intensity boundary

for i=1:Ncells
    DumKdi=Kdi(:,i)-LB;
    DumKdi(DumKdi<0)=0;
    DumKdi(DumKdi>0)=DumKdi(DumKdi>0)+LB;
    FKdi(:,i)=DumKdi;
    
    DumKdiR1=KdiR1(:,i)-LB;
    DumKdiR1(DumKdiR1<0)=0;
    DumKdiR1(DumKdiR1>0)=DumKdiR1(DumKdiR1>0)+LB;
    FKdiR1(:,i)=DumKdiR1;
    
    DumKdiR2=KdiR2(:,i)-LB;
    DumKdiR2(DumKdiR2<0)=0;
    DumKdiR2(DumKdiR2>0)=DumKdiR2(DumKdiR2>0)+LB;
    FKdiR2(:,i)=DumKdiR2;
    
    % Filtering -- inter spot distances
    
    for j=1:length(Kdx(:,i))
    Ddiff=abs(Kdx(j,i)-KdxR1(j,i));
    Ddiff2=abs(Kdx(j,i)-KdxR2(j,i));
    if Ddiff<Dcrit
        FKdi(j,i)=FKdi(j,i)+FKdiR1(j,i);
        FKdiR1(j,i)=0;
        SN(j,i)=1;
    else
        FKdi(j,i)=FKdi(j,i)+FKdiR1(j,i);
        FKdiR1(j,i)=0;
        SN(j,i)=2;
    end
    
    if Ddiff2<Dcrit
        FKdi(j,i)=FKdi(j,i)+FKdiR2(j,i);
        FKdiR2(j,i)=0;

    else
        FKdi(j,i)=FKdi(j,i)+FKdiR2(j,i);
        FKdiR2(j,i)=0;

    end
        
    end
end

%% Plot Spotnumber distribution

load('SNumber.mat')

X4=(1:length(SNumber))/length(SNumber);

plot(X4,SNumber,'or')
axis([0 1 0.5 2.5])
title('Spot Number vs. Time','FontSize',20)
xlabel('Normalized Cell Time(-)','FontSize',18)
ylabel('Spot Number (-)','FontSize',18)

%% Test plots
% plots 10 8 6 4 3 1
% 
for i=3
figure(i)
hold on
 plot(xccd,FKdi(:,i),'ob','LineWidth',4);
 plot(xccd,FKdiR1(:,i),'or','LineWidth',4);
% plot(xccd,FKdiR2(:,i),'oc','LineWidth',4);
%plot(xccd,FKdi(:,i)+FKdiR1(:,i)+FKdiR2(:,i),'ok','LineWidth',4);
hold off

% figure(i+1)
% hold on 
% plot(xccd,abs(Kdx(:,i)-KdxR2(:,i)),'b','LineWidth',2)
% hold off
end

%% Mean CALC
MeanMat=[];
MeanMatR1=[];
MeanMatR2=[];

for i=[1 3 4 6 8 10]
    MK(:,i)=FKdi(:,i);
    MKR1(:,i)=FKdiR1(:,i);
    MKR2(:,i)=FKdiR2(:,i);
    
    MeanMat=[MeanMat MK(:,i)];
    MeanMatR1=[MeanMatR1 MKR1(:,i)];
    MeanMatR2=[MeanMatR2 MKR2(:,i)];
end
%% plot Mean

MEAN=mean(MeanMat,2);
STDMEAN=std(MEAN);
MEANR1=mean(MeanMatR1,2);
MEANR2=mean(MeanMatR2,2);

load('MEanEdit.mat');

X=(1:length(MEANedit(:,2)))/length(MEANedit(:,2));
X2=(1:length(MEAN))/length(MEAN);

STDMEANedit=std(MEANedit(:,2));

hold on
plot(X,MEANedit(:,2),'-ob');
plot(X,MEANedit(:,2)+STDMEANedit,':b');
plot(X,MEANedit(:,2)-STDMEANedit,':b');
% plot(X,MEANR1,'-or');
hold off 
title('Mean Integrated Intensity vs. Time','FontSize',20)
xlabel('Normalized Cell Time (-)','FontSize',20,'FontWeight','bold')
ylabel('Integrated Intensity (-)','FontSize',20,'FontWeight','bold')


%% 
i=9;
hold on
% plot(X,Ki(:,1)+KiR1(:,1));
% plot(X,Ki(:,2)+KiR1(:,2));
% plot(X,Ki(:,3));
%   plot(X,Ki(:,i),'r');
  plot(X2,Ki(:,i)+KiR1(:,i),'b');
%   plot(X,Ki(:,i)+KiR1(:,i)+KiR2(:,i),'k');
% plot(X,Ki(:,6));
% plot(X,Ki(:,7));
hold off

 %% DnaN
%  M(40:45,1)=M(40:45,1)/1.1;
%  M(45:50,1)=M(45:50,1)/1.2;
%  M(50:55,1)=M(50:55,1)/1.2;
%  M(55:60,1)=M(55:60,1)/1.2;
 
hold on
plot(X2,M(:,1))
plot(X2,(M(:,1)+M(:,6)/1.5))
plot(X2,(M(:,1)-M(:,6)/1.5))
hold off