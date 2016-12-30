clear all
close
clc

folder='/Users/rleeuw/Work/Data/160205_Agar_Data';
slash = '/';
exps=[1 2  3 4 5 7 8 9];
Intensityval = [700, 300, 700]; %[CFP YFP RFP]
umperpx=0.159;

Lcfp=[];    Lyfp=[];    Lrfp=[];
Pcfp=[];    Pyfp=[];    Prfp=[];
Icfp=[];    Iyfp=[];    Irfp=[];
Fcfp=[];    Fyfp=[];    Frfp=[];
ncfp=[];    nyfp=[];    nrfp=[];
celllength = [];

yfp.filterval=Intensityval(2)*2;
cfp.filterval=Intensityval(1)*2;
rfp.filterval=Intensityval(3)*2;
    
j=1; 
for i=exps;
    E{j}=load(strcat(folder,slash,num2str(i),slash,'Results.mat')); 
    j=j+1;
end

allCFP_L = [];

Nexp=size(E,2);


for i=1:Nexp
    
    Ncells{i}=size(E{i}.DataStruct,2);
    
    
    for j=1:Ncells{i} 
        
        if ~isempty(E{i}.DataStruct(1,j).Lnorm)        
            LNormCFP{i,j}=E{i}.DataStruct(1,j).Lnorm;
        else
            LNormCFP{i,j}=0;
        end
        
        CellLength{i,j}=E{i}.DataStruct(1,j).CellLength;
        
        
        CFPld{i,j}=E{i}.DataStruct(1,j).ld;
        YFPld{i,j}=E{i}.DataStruct(2,j).ld;
        RFPld{i,j}=E{i}.DataStruct(3,j).ld;
                
        NspotsCFP=size(CFPld{i,j},2);
        NspotsYFP=size(YFPld{i,j},2);
        NspotsRFP=size(RFPld{i,j},2);        
     
        if NspotsCFP==0
            CFPld{i,j}{1}=[];
        else
            for k=1:NspotsCFP
                Lcfp=[Lcfp CellLength{i,j}];
                Pcfp=[Pcfp CFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Icfp=[Icfp CFPld{i,j}{k}(1,1)/Intensityval(1)];
                Fcfp=[Fcfp CFPld{i,j}{k}(1,7)];
            end
        end
        ncfp = [ncfp NspotsCFP];
        celllength = [celllength CellLength{i,j}];
        
        if NspotsYFP==0
            YFPld{i,j}{1}=[];
        else
            for k=1:NspotsYFP
                Lyfp=[Lyfp CellLength{i,j}];
                Pyfp=[Pyfp YFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Iyfp=[Iyfp YFPld{i,j}{k}(1,1)/Intensityval(2)];
                Fyfp=[Fyfp YFPld{i,j}{k}(1,7)];
            end
        end
        nyfp = [nyfp NspotsYFP];
        
        if NspotsRFP==0
            RFPld{i,j}{1}=[];
        else
            for k=1:NspotsRFP            
                Lrfp=[Lrfp CellLength{i,j}];
                Prfp=[Prfp RFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Irfp=[Irfp RFPld{i,j}{k}(1,1)/Intensityval(3)];
                Frfp=[Frfp RFPld{i,j}{k}(1,7)];
            end
        end
        nrfp = [nrfp NspotsRFP];
    end
end

% tozero = Iyfp<yfp.filterval;
% Iyfp(tozero) = 0;
% Lyfp(tozero) = 0;
% Pyfp(tozero) = 0;
% Fyfp(tozero) = 0;


%% Position vs. cell length
% CFP

fig1 = figure(1);
set(fig1,'Position',[20,300,1800,500])


subplot(1,3,1)
hold on
scatter(single(Lcfp),Pcfp,Icfp,'b','filled');
myfit=polyfit(Lcfp,Pcfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
% plot(x,y,'r','LineWidth',5)
xlabel('Cell Length'); ylabel('Normalized osition of spot in cell'); 
title('Agar data: CFP')
hold off
axis([12 43 -0.1 1.1])
set(gca,'FontSize',16)


% YFP

subplot(1,3,2);
hold on
scatter(single(Lyfp),Pyfp,Iyfp,'m','filled');
myfit=polyfit(Lyfp,Pyfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
xlabel('Cell Length'); ylabel('Normalized osition of spot in cell'); 
title('Agar data: YFP')
hold off
axis([12 43 -0.1 1.1])
set(gca,'FontSize',16)

% RFP

subplot(1,3,3)
hold on
scatter(single(Lrfp),Prfp,Irfp,'r','filled');
myfit=polyfit(Lrfp,Prfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
xlabel('Cell Length'); ylabel('Normalized osition of spot in cell'); 
title('Agar data: RFP')
% plot(x,y,'k','LineWidth',5)
hold off
axis([12 43 -0.1 1.1])
set(gca,'FontSize',16)

%% Intensity vs. position

% CFP

fig2 = figure(2);
set(fig2,'Position',[20,300,1800,500])
subplot(1,3,1)

hold on
scatter(Pcfp,Icfp,'b','x');
myfit=polyfit(Pcfp,Icfp,4);
x=0:0.001:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Agar data: CFP')
hold off
axis([0 1 -0.1 90])
set(gca,'FontSize',16)


% YFP

subplot(1,3,2)
hold on
scatter(Pyfp,Iyfp,'m','x');
myfit=polyfit(Pyfp,Iyfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Agar data: YFP')
hold off
axis([0 1 -0.1 90])
set(gca,'FontSize',16)


% RFP

subplot(1,3,3)
hold on
scatter(Prfp,Irfp,'r','x');
myfit=polyfit(Prfp,Irfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Agar data: RFP')
hold off
axis([0 1 -0.1 90])
set(gca,'FontSize',16)


%% Numspots vs. position

bins = 15;
thisedge = (0:bins)/bins;

fig3 = figure(3);
set(fig3,'Position',[20,300,1800,500])

% CFP

subplot(1,3,1)
[numbin,edges] = histcounts(Pcfp,thisedge);
norm = max(numbin)/80;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Pcfp,Icfp,'b','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Agar data: CFP')
axis([0 1 -0.1 90])
set(gca,'FontSize',16)

% YFP

subplot(1,3,2)
[numbin,edges] = histcounts(Pyfp,thisedge);
norm = max(numbin)/80;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Pyfp,Iyfp,'m','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Agar data: YFP')
axis([0 1 -0.1 90])
set(gca,'FontSize',16)

% RFP

subplot(1,3,3)
[numbin,edges] = histcounts(Prfp,thisedge);
norm = max(numbin)/80;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Prfp,Irfp,'r','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Agar data: RFP')
axis([0 1 -0.1 90])
set(gca,'FontSize',16)


%% Numspots vs. position & cell length

bins = 15;
thisedge2{1} = linspace(min(Lcfp),max(Lcfp),bins+1);
thisedge2{2} = (0:bins)/bins;

fig4 = figure(4);
set(fig4,'Position',[20,300,1800,500])
fig5 = figure(5);
set(fig5,'Position',[20,300,1800,500])

% CFP

subplot(1,3,1)
Numcfp(1,:) = Lcfp;
Numcfp(2,:) = Pcfp;

%Filter on Intensity
ICFP=Icfp*Intensityval(1);
[bin_cfp,idx_cfp]=find(ICFP>cfp.filterval); 

NumCFP(:,idx_cfp)=Numcfp(:,idx_cfp);
NumCFPnz(1,:)=nonzeros(NumCFP(1,:));
NumCFPnz(2,:)=nonzeros(NumCFP(2,:));

FilteredSignalCFP=size(NumCFPnz,2)/size(ICFP,2); % Percentage of total signal

figure(4)
subplot(1,3,1)
Heatmap = hist3(NumCFPnz','Edges',thisedge2);
pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
colormap(fig4,jet) % heat map
xlabel('Cell Length'); ylabel('Position in Cell');
title('Agar data: CFP');
grid on
set(gca,'FontSize',16)

figure(5)
subplot(1,3,1)
hold on
hist3(NumCFPnz','Edges',thisedge2)
colormap(fig5,jet) % heat map
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Agar data: CFP');
grid off
hold off
axis([min(Lcfp),max(Lcfp),0,1])
view(3)
set(gca,'FontSize',16)

% YFP
Numyfp(1,:) = Lyfp;
Numyfp(2,:) = Pyfp;

%Filter on Intensity
IYFP=Iyfp*Intensityval(2);
[bin_yfp,idx_yfp]=find(IYFP>yfp.filterval); 

NumYFP(:,idx_yfp)=Numyfp(:,idx_yfp);
NumYFPnz(1,:)=nonzeros(NumYFP(1,:));
NumYFPnz(2,:)=nonzeros(NumYFP(2,:));

FilteredSignalYFP=size(NumYFPnz,2)/size(IYFP,2); % Percentage of total signal

figure(4)
subplot(1,3,2)
Heatmap = hist3(NumYFPnz','Edges',thisedge2);
pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
xlabel('Cell Length'); ylabel('Position in Cell');
title('Agar data: YFP');
grid on
set(gca,'FontSize',16)

figure(5)
subplot(1,3,2)
hold on
hist3(NumYFPnz','Edges',thisedge2)
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Agar data: YFP');
grid off
hold off
axis([min(Lyfp),max(Lyfp),0,1])
view(3)
set(gca,'FontSize',16)

% RFP
Numrfp(1,:) = Lrfp;
Numrfp(2,:) = Prfp;

%Filter on Intensity
IRFP=Irfp*Intensityval(3);
[bin_rfp,idx_rfp]=find(IRFP>rfp.filterval); 

NumRFP(:,idx_rfp)=Numrfp(:,idx_rfp);

NumRFPnz(1,:)=nonzeros(NumRFP(1,:));
NumRFPnz(2,:)=nonzeros(NumRFP(2,:));

FilteredSignalRFP=size(NumRFPnz,2)/size(IRFP,2);

figure(4)
subplot(1,3,3)
Heatmap = hist3(NumRFPnz','Edges',thisedge2);
h = pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
xlabel('Cell Length'); ylabel('Position in Cell');
title('Agar data: RFP');
grid on
set(gca,'FontSize',16)

figure(5)
subplot(1,3,3)
hold on
hist3(NumRFPnz','Edges',thisedge2)
colormap(jet) % heat map
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Agar data: RFP');
grid off
hold off
axis([min(Lrfp),max(Lrfp),0,1])
view(3)
set(gca,'FontSize',16)

%% Full cell intensity vs. celllength

clear plotcfp plotyfp plotrfp
fig6 = figure(6);
set(fig6,'Position',[20,300,1800,500])

% CFP
bins1 = 7;
bins2 = 50;
thisedge3{1} = linspace(min(Lcfp),max(Lcfp),bins1+1);
thisedge3{2} = linspace(min(Fcfp),max(Fcfp),bins2+1);

plotcfp(1,:) = Lcfp; %length
plotcfp(2,:) = Fcfp; %full cell intensity
plotcfp(3,:) = Icfp*Intensityval(1); %spot intensity

plotcfp = unique(plotcfp','rows')';
% 
% [N_full,Edges_full,mid_full,loc_full]=histcn([plotcfp(1,:)' plotcfp(2,:)'],thisedge3{1},thisedge3{2});
% [N_spot,Edges_spot,mid_spot,loc_spot]=histcn([plotcfp(1,:)' plotcfp(3,:)'],thisedge3{1},thisedge3{2});
% 
% for i=1:size(N_full,1);
% CFP_binned(i,:)=N_full(i,:).*Edges_full{2};
% Y=CFP_binned(i,:)';
% Y=double(Y);
% X=1:length(CFP_binned(i,:));
% display(strcat({'Gaussian fit of CFP column '},num2str(i),{' of '},num2str(size(N_full,1))));
% f{i}=fit(X',Y,'gauss1');
% peak(i)=f{i}.b1;
% sigma(i)=f{i}.c1;
% peakfloored(i)=floor(peak(i));
% peakceiled(i)=ceil(peak(i));
% IntVal(i)=((Edges_full{2}(peakceiled(i))-Edges_full{2}(peakfloored(i))))*(peakceiled(i)-peak(i))+Edges_full{2}(peakfloored(i)); %slope * peakposition, because linear.
% end


subplot(1,3,1)
hold on
scatter(plotcfp(1,:),plotcfp(2,:),'b','o','filled');
scatter(plotcfp(1,:),plotcfp(3,:),'r','o','filled');
myfit=polyfit(plotcfp(1,:),plotcfp(2,:),1);
myfit2=polyfit(plotcfp(1,:),plotcfp(3,:),1);
x=12:0.1:43;
y=polyval(myfit,x);
y2=polyval(myfit2,x);
plot(x,y,'b','LineWidth',3)
plot(x,y2,'r','LineWidth',3)
xlabel('Cell Length (px)'); ylabel('Intensity (-)'); 
title('TetR stoichiometry vs. Length')
hold off
axis([12 43 -0.1 2*10^5])
set(gca,'FontSize',16)
legend('Full cell intensity','Spot intensity')

clear N_full
% YFP


plotyfp(1,:) = Lyfp;
plotyfp(2,:) = Fyfp;
plotyfp(3,:) = Iyfp*Intensityval(2);

plotyfp = unique(plotyfp','rows')';

[N_full,Edges_full,mid_full,loc_full]=histcn([plotyfp(1,:)' plotyfp(2,:)'],thisedge3{1},thisedge3{2});
[N_spot,Edges_spot,mid_spot,loc_spot]=histcn([plotyfp(1,:)' plotyfp(3,:)'],thisedge3{1},thisedge3{2});

for i=1:size(N_full,1);
YFP_binned(i,:)=N_full(i,:).*Edges_full{2};
Y=YFP_binned(i,:)';
Y=double(Y);
X=1:length(YFP_binned(i,:));
display(strcat({'Gaussian fit of CFP column '},num2str(i),{' of '},num2str(size(N_full,1))));
f{i}=fit(X',Y,'gauss1');
peak(i)=f{i}.b1;
sigma(i)=f{i}.c1;
peakfloored(i)=floor(peak(i));
peakceiled(i)=ceil(peak(i));
IntVal(i)=((Edges_full{2}(peakceiled(i))-Edges_full{2}(peakfloored(i))))*(peakceiled(i)-peak(i))+Edges_full{2}(peakfloored(i)); %slope * peakposition, because linear.
end
plot(Edges_full{1},IntVal);

subplot(1,3,2)
hold on
scatter(plotyfp(1,:),plotyfp(2,:),'b','o','filled');
scatter(plotyfp(1,:),plotyfp(3,:),'r','o','filled');
myfit=polyfit(plotyfp(1,:),plotyfp(2,:),1);
myfit2=polyfit(plotyfp(1,:),plotyfp(3,:),1);
x=12:0.1:43;
y=polyval(myfit,x);
y2=polyval(myfit2,x);
plot(x,y,'b','LineWidth',3)
% plot(x,y2,'r','LineWidth',3)
xlabel('Cell Length'); ylabel('Intensity'); 
title('Tus Stoichiometry vs. Length')
hold off
axis([12 43 -0.1 2*10^5])
set(gca,'FontSize',16)


% RFP

plotrfp(1,:) = Lrfp;
plotrfp(2,:) = Frfp;
plotrfp(3,:) = Irfp*Intensityval(3);
plotrfp = unique(plotrfp','rows')';

subplot(1,3,3)
hold on
scatter(plotrfp(1,:),plotrfp(2,:),'b','o','filled');
scatter(plotrfp(1,:),plotrfp(3,:),'r','o','filled');
myfit=polyfit(plotrfp(1,:),plotrfp(2,:),1);
myfit2=polyfit(plotrfp(1,:),plotrfp(3,:),1);
x=12:0.1:43;
y=polyval(myfit,x);
y2=polyval(myfit2,x);
plot(x,y,'b','LineWidth',3)
plot(x,y2,'r','LineWidth',3)
xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
title('DnaN Stoichiometry vs. Length')
hold off
axis([12 43 -0.1 4*10^5])
set(gca,'FontSize',16)

%% Numspots/cell vs. cell length
fig1 = figure(7);

for i=1:max(nrfp);
p{i} = find(nrfp == i);
m(i) = mean(celllength(p{i}));
Spotstd(i)=std(celllength(p{i}));
end

myfit=polyfit(m,1:4,1);
x=[12,43];
y=polyval(myfit,x);

Y=linspace(1,max(nrfp),max(nrfp));

Xl=m-Spotstd;
Xr=m+Spotstd;

hold on
scatter(celllength,nrfp,'m')
plot(m,Y,'b--','LineWidth',3)
plot(Xl,Y,'b','LineWidth',1)
plot(Xr,Y,'b','LineWidth',1)
axis([12, 43, 0, 5])
xlabel('Cell length'); ylabel('Number of spots per cell')
title('Numspots/cell vs. cell length for YFP')
hold off
set(gca,'FontSize',16)
legend('Spot vs Length points','Mean value','std')