clear 
close all
clc

folder='C:\Users\water\Documents\GitHub\Data\Target Data\Kymo Data\Results';
slash = '\';
channels = [1,2];
Intensityval = [350,60,200];

Lcfp=[];    Lyfp=[];    Lrfp=[];
Pcfp=[];    Pyfp=[];    Prfp=[];
Icfp=[];    Iyfp=[];    Irfp=[];
Fcfp=[];    Fyfp=[];    Frfp=[];

for i=channels;
    E{i}=load(strcat(folder,slash,'Results_Ch',num2str(i),'.mat')); 
end

allCFP_L = [];


for i=channels
    
    Ncells{i}=size(E{i}.DataStruct,2);
    
    
    for j=1:Ncells{i}     
        
        CFPx{i}{j}=E{i}.DataStruct(1,j).x;
        YFPx{i}{j}=E{i}.DataStruct(2,j).x;
        RFPx{i}{j}=E{i}.DataStruct(3,j).x;
                
        NspotsCFP=size(CFPx{i}{j},2);
        NspotsYFP=size(YFPx{i}{j},2);
        NspotsRFP=size(RFPx{i}{j},2);
        
        frames{i}{j} = size(CFPx{i}{j}{1},1);
     
        if NspotsCFP==0
            CFPx{i}{j}{1}=[];
        else
            for k=1:NspotsCFP
                for h = 1:frames{i}{j};
                    length = size(E{i}.DataStruct(1,j).ydatacrpdR1{h,k},2);
                    Lcfp=[Lcfp h/frames{i}{j}];
                    Pcfp=[Pcfp CFPx{i}{j}{k}(h,2)/length];
                    Icfp=[Icfp CFPx{i}{j}{k}(h,1)/Intensityval(1)];
                    Fcfp=[Fcfp CFPx{i}{j}{k}(h,7)];
                end
            end
        end
        
        if NspotsYFP==0
            YFPx{i}{j}{1}=[];
        else
            for k=1:NspotsYFP
                for h = 1:frames{i}{j};
                    length = size(E{i}.DataStruct(2,j).ydatacrpdR1{h,k},2);
                    Lyfp=[Lyfp h/frames{i}{j}];
                    Pyfp=[Pyfp YFPx{i}{j}{k}(h,2)/length];
                    Iyfp=[Iyfp YFPx{i}{j}{k}(h,1)/Intensityval(2)];
                    Fyfp=[Fyfp YFPx{i}{j}{k}(h,7)];
                end
            end
        end
        
        if NspotsRFP==0
            RFPx{i}{j}{1}=[];
        else
            for k=1:NspotsRFP
                for h = 1:frames{i}{j};
                    length = size(E{i}.DataStruct(3,j).ydatacrpdR1{h,k},2);
                    Lrfp=[Lrfp h/frames{i}{j}];
                    Prfp=[Prfp RFPx{i}{j}{k}(h,2)/length];
                    Irfp=[Irfp RFPx{i}{j}{k}(h,1)/Intensityval(3)];
                    Frfp=[Frfp RFPx{i}{j}{k}(h,7)];
                end
            end
        end
    end
end

Yremove = Iyfp == 0;
Iyfp(Yremove) = [];
Lyfp(Yremove) = [];
Pyfp(Yremove) = [];
Fyfp(Yremove) = [];

Rremove = Irfp == 0;
Irfp(Rremove) = [];
Lrfp(Rremove) = [];
Prfp(Rremove) = [];
Frfp(Rremove) = [];

Cremove = Icfp == 0;
Icfp(Cremove) = [];
Lcfp(Cremove) = [];
Pcfp(Cremove) = [];
Fcfp(Cremove) = [];

%% Position vs. cell length
% CFP

fig1 = figure(1);
set(fig1,'Position',[20,300,1800,500])
subplot(1,3,1)
hold on
scatter(single(Lcfp),Pcfp,Icfp,'b','filled');
% myfit=polyfit(Acfp,Bcfp,4);
% x=0:0.01:1;
% y=polyval(myfit,x);
% plot(x,y,'r','LineWidth',5)
hold off
xlabel('Cell Length'); ylabel('Normalized position of spot in cell'); 
title('Kymo data: CFP')
axis([0 1 0 1])

% YFP

figure(1)
subplot(1,3,2)
hold on
scatter(single(Lyfp),Pyfp,Iyfp,'m','filled');
% myfit=polyfit(Ayfp,Byfp,4);
% x=0:0.01:1;
% y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
hold off
xlabel('Cell Length'); ylabel('Normalized position of spot in cell'); 
title('Kymo data: YFP')
axis([0 1 0 1])

% RFP

figure(1)
subplot(1,3,3)
hold on
scatter(single(Lrfp),Prfp,Irfp,'r','filled');
% myfit=polyfit(Arfp,Brfp,4);
% x=0:0.01:1;
% y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
hold off
xlabel('Cell Length'); ylabel('Normalized position of spot in cell'); 
title('Kymo data: RFP')
axis([0 1 0 1])
        

%% Intensity vs. position

% CFP

fig2 = figure(2);
set(fig2,'Position',[20,300,1800,500])
subplot(1,3,1)
hold on
scatter(Pcfp,Icfp,'b','x');
myfit=polyfit(Pcfp,Icfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Kymo data: CFP')
hold off
axis([0 1 -0.1 35])

% YFP

figure(2)
subplot(1,3,2)
hold on
scatter(Pyfp,Iyfp,'m','x');
myfit=polyfit(Pyfp,Iyfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Kymo data: YFP')
hold off
axis([0 1 -0.1 35])

% RFP

figure(2)
subplot(1,3,3)
hold on
scatter(Prfp,Irfp,'r','x');
myfit=polyfit(Prfp,Irfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Kymo data: RFP')
hold off
axis([0 1 -0.1 35])

%% Numspots vs. position
% CFP
fig3 = figure(3);
set(fig3,'Position',[20,300,1800,500])
subplot(1,3,1)

[numbin,edges] = histcounts(Pcfp,20);
norm = max(numbin)/35;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Pcfp,Icfp,'b','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Kymo data: CFP')
axis([0 1 -0.1 40])

% YFP
figure(3);
subplot(1,3,2)

[numbin,edges] = histcounts(Pyfp,20);
norm = max(numbin)/35;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Pyfp,Iyfp,'m','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Kymo data: YFP')
axis([0 1 -0.1 40])

% RFP
figure(3);
subplot(1,3,3)

[numbin,edges] = histcounts(Prfp,20);
norm = max(numbin)/35;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Prfp,Irfp,'r','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)'); 
title('Kymo data: RFP')
axis([0 1 -0.1 40])

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

figure(4)
subplot(1,3,1)
Heatmap = hist3(Numcfp','Edges',thisedge2);
pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
colormap(hot) % heat map
xlabel('Cell Length'); ylabel('Position in Cell');
title('Kymo data: CFP');
grid on

figure(5)
subplot(1,3,1)
hold on
hist3(Numcfp','Edges',thisedge2)
colormap(hot) % heat map
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Kymo data: CFP');
grid off
hold off
axis([min(Lcfp),max(Lcfp),0,1])
view(3)

% YFP
Numyfp(1,:) = Lyfp;
Numyfp(2,:) = Pyfp;

figure(4)
subplot(1,3,2)
Heatmap = hist3(Numyfp','Edges',thisedge2);
pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
colormap(hot) % heat map
xlabel('Cell Length'); ylabel('Position in Cell');
title('Kymo data: YFP');
grid on

figure(5)
subplot(1,3,2)
hold on
hist3(Numyfp','Edges',thisedge2)
colormap(hot) % heat map
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Kymo data: YFP');
grid off
hold off
axis([min(Lyfp),max(Lyfp),0,1])
view(3)

% RFP
Numrfp(1,:) = Lrfp;
Numrfp(2,:) = Prfp;

figure(4)
subplot(1,3,3)
Heatmap = hist3(Numrfp','Edges',thisedge2);
h = pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
colormap(hot) % heat map
xlabel('Cell Length'); ylabel('Position in Cell');
title('Kymo data: RFP');
grid on

figure(5)
subplot(1,3,3)
hold on
hist3(Numrfp','Edges',thisedge2)
colormap(hot) % heat map
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('Cell Length'); ylabel('Position in Cell');zlabel('Amount of spots')
title('Kymo data: RFP');
grid off
hold off
axis([min(Lrfp),max(Lrfp),0,1])
view(3)


%% Full cell intensity vs. celllength

clear plotcfp plotyfp plotrfp
fig6 = figure(6);
set(fig6,'Position',[20,300,1800,500])

% YFP

plotcfp(1,:) = Lcfp;
plotcfp(2,:) = Fcfp/max(Fcfp);
plotcfp = unique(plotcfp','rows')';

subplot(1,3,1)
hold on
scatter(plotcfp(1,:),plotcfp(2,:),'b','x');
myfit=polyfit(plotcfp(1,:),plotcfp(2,:),1);
x=0:0.01:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
title('Kymo data: CFP')
hold off
axis([0 1 -0.1 1.1])

% YFP

plotyfp(1,:) = Lyfp;
plotyfp(2,:) = Fyfp/max(Fyfp);
plotyfp = unique(plotyfp','rows')';

subplot(1,3,2)
hold on
scatter(plotyfp(1,:),plotyfp(2,:),'m','x');
myfit=polyfit(plotyfp(1,:),plotyfp(2,:),1);
x=0:0.01:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
title('Kymo data: YFP')
hold off
axis([0 1 -0.1 1.1])


% RFP

plotrfp(1,:) = Lrfp;
plotrfp(2,:) = Frfp/max(Frfp);
plotrfp = unique(plotrfp','rows')';

subplot(1,3,3)
hold on
scatter(plotrfp(1,:),plotrfp(2,:),'r','x');
myfit=polyfit(plotrfp(1,:),plotrfp(2,:),1);
x=0:0.01:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
title('Kymo data: RFP')
hold off
axis([0 1 -0.1 1.1])