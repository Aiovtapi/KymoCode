GaussCalcs

%%
LBd=5000;
LB=800;

FKdi=zeros(MaxBacLifed,Ncells);
FKdiR1=zeros(MaxBacLifed,Ncells);
FKdiR2=zeros(MaxBacLifed,Ncells);

FKi=zeros(MaxBacLife,Ncells);
FKiR1=zeros(MaxBacLife,Ncells);
FKiR2=zeros(MaxBacLife,Ncells);

MK=zeros(MaxBacLifed,Ncells);

Dcrit=0; % minimal distance between spots for identification

% Filtering -- lower intensity boundary

for i=1:Ncells
    
    %   dif ----------------
    DumKdi=Kdi(:,i)-LBd;
    DumKdi(DumKdi<0)=0;
    DumKdi(DumKdi>0)=DumKdi(DumKdi>0)+LBd;
    FKdi(:,i)=DumKdi;
    
    DumKdiR1=KdiR1(:,i)-LBd;
    DumKdiR1(DumKdiR1<0)=0;
    DumKdiR1(DumKdiR1>0)=DumKdiR1(DumKdiR1>0)+LBd;
    FKdiR1(:,i)=DumKdiR1;
    
    DumKdiR2=KdiR2(:,i)-LBd;
    DumKdiR2(DumKdiR2<0)=0;
    DumKdiR2(DumKdiR2>0)=DumKdiR2(DumKdiR2>0)+LBd;
    FKdiR2(:,i)=DumKdiR2;
    
    % Filtering -- interspot distances
    for j=1:length(Kdx(:,i))
    Ddiff=abs(Kdx(j,i)-KdxR1(j,i));
    Ddiff2=abs(Kdx(j,i)-KdxR2(j,i));
    Ddiff3=abs(KdxR1(j,i)-KdxR2(j,i));
    if Ddiff<Dcrit
        FKdi(j,i)=FKdi(j,i)+FKdiR1(j,i);
        FKdiR1(j,i)=0;
    end
    if Ddiff2<Dcrit
        FKdi(j,i)=FKdi(j,i)+FKdiR2(j,i);
        FKdiR2(j,i)=0;
    end
    if Ddiff3<Dcrit && Ddiff2>Dcrit
        FKdiR1(j,i)=FKdiR1(j,i)+FKdiR2(j,i);
        FKdiR2(j,i)=0;
    end
    end
    
    % -------------- TUS --------------
    
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
    
    % Filtering -- interspot distances
    for j=1:length(Kx(:,i))
    Ddiff=abs(Kx(j,i)-KxR1(j,i));
    Ddiff2=abs(Kx(j,i)-KxR2(j,i));
    Ddiff3=abs(KxR1(j,i)-KxR2(j,i));
    if Ddiff<Dcrit
        FKi(j,i)=FKi(j,i)+FKiR1(j,i);
        FKiR1(j,i)=0;
    end
    if Ddiff2<Dcrit
        FKi(j,i)=FKi(j,i)+FKiR2(j,i);
        FKiR2(j,i)=0;
    end
    if Ddiff3<Dcrit && Ddiff2>Dcrit
        FKiR1(j,i)=FKiR1(j,i)+FKiR2(j,i);
        FKiR2(j,i)=0;
    end
    end
end
X=(1:length(Ki(:,1)))/length(Ki(:,1));
Xd=(1:length(Kdi(:,1)))/length(Kdi(:,1));

Kdxleft=[];
Kdxright=[];
KdxleftR1=[];
KdxrightR1=[];
Kxleft=[];
Kxright=[];
KxleftR1=[];
KxrightR1=[];
KxleftR2=[];
KxrightR2=[];
KxleftR3=[];
KxrightR3=[];

for i=1:Ncells
if Kdx(end,i)<=0.5
    Kdxleft=[Kdxleft Kdx(:,i)];
else
    Kdxright=[Kdxright Kdx(:,i)];
end

if KdxR1(end,i)<=0.5
    KdxleftR1=[KdxleftR1 KdxR1(:,i)];
else
    KdxrightR1=[KdxrightR1 KdxR1(:,i)];
end

if Kx(end,i)<=0.5 
    Kxleft=[Kxleft Kx(:,i)];
else
    Kxright=[Kxright Kx(:,i)];
end

if KxR1(end,i)<=0.5 
    KxleftR1=[KxleftR1 KxR1(:,i)];
else
    KxrightR1=[KxrightR1 KxR1(:,i)];
end

if KxR2(end,i)<=0.5 
    KxleftR2=[KxleftR2 KxR2(:,i)];
else
    KxrightR2=[KxrightR2 KxR2(:,i)];;
end

if KxR3(end,i)<=0.5 
    KxleftR3=[KxleftR3 KxR3(:,i)];
else
    KxrightR3=[KxrightR3 KxR3(:,i)];
end
end


%% Some more calcs
MeanKdxleft=mean(Kdxleft,2);
MeanKdxright=mean(Kdxright,2);

MeanKxleft=mean(Kxleft,2);
MeanKxright=mean(Kxright,2);

MeanKxleftR1=mean(KxleftR1,2);
MeanKxrightR1=mean(KxrightR1,2);

MeanKxleftR2=mean(KxleftR2,2);
MeanKxrightR2=mean(KxrightR2,2);

%% Plots

Fac=1/(length(Xd)/(MeanCellLifed+9));

hold on
plot(Md(1:47,2),Xd(1:47)/Fac,'r')
plot(MeanKdxleft(47:94),Xd(47:94)/Fac,'r')
plot(MeanKdxright(47:94),Xd(47:94)/Fac,'r')
plot(MR2(1:15,2),X(1:15)/Fac,'b')
plot(MeanKxrightR2(15:30),X(15:30)/Fac,'b');
plot(mean([MeanKxleftR2(15:30) MeanKxleftR1(15:30) MeanKxleft(15:30)],2),X(15:30)/Fac,'b');
hold off
line([0 1],[1 1],'LineWidth',2)
line([0.5 0.5],[1 2])
axis([0 1 0 2])
title('Tus and dif localization w.r.t. cell time')
xlabel('Normalized Cell Position')
ylabel('Normalized Cell Time')

%% Tus w.r.t. dif

Fac=1/(length(Xd)/(MeanCellLifed+10));

figure(8)
hold on
line([0 1],[1 1],'LineWidth',2)
line([0.5 0.5],[1 2])
% hline2=refline([0 MeanCellLifeT*5]);
% hline2.Color='k';
axis([0 1 0 2])
hold off