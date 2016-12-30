clear 
close all
clc

folder='/Users/rleeuw/Work/Data/141230_Kymo_Results';
slash = '/';
channels = [1,2];
Intensityval = [350,60,200];

Lcfp=[];    Lyfp=[];    Lrfp=[];
Pcfp=[];    Pyfp=[];    Prfp=[];
Icfp=[];    Iyfp=[];    Irfp=[];
Fcfp=[];    Fyfp=[];    Frfp=[];
ncfp=[];    nyfp=[];    nrfp=[];

YFP.filterval=700;
RFP.filterval=700;
CFP.filterval=1000;



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
                    IcfpM{i}{j}(h,k)=CFPx{i}{j}{k}(h,1);
                    IcfpM{i}{j}(IcfpM{i}{j}<CFP.filterval)=0;
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
                    lyfp{i}{j}(h)=h/frames{i}{j};
                    Pyfp=[Pyfp YFPx{i}{j}{k}(h,2)/length];
                    Iyfp=[Iyfp YFPx{i}{j}{k}(h,1)/Intensityval(2)];
                    IyfpM{i}{j}(h,k)=YFPx{i}{j}{k}(h,1);
                    IyfpM{i}{j}(IyfpM{i}{j}<YFP.filterval)=0;
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
                    lrfp{i}{j}=h/frames{i}{j};
                    Prfp=[Prfp RFPx{i}{j}{k}(h,2)/length];
                    Irfp=[Irfp RFPx{i}{j}{k}(h,1)/Intensityval(3)];                   
                    IrfpM{i}{j}(h,k)=RFPx{i}{j}{k}(h,1);
                    IrfpM{i}{j}(IrfpM{i}{j}<RFP.filterval)=0;
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

IYFP=Iyfp*Intensityval(2);
ICFP=Icfp*Intensityval(1);
IRFP=Irfp*Intensityval(3);

%Total spot intensities
%Spot number

for i=channels
    for j=1:Ncells{i}
        
        IyfpSpotTotal{i,j}=sum(IyfpM{i}{j},2);
        IcfpSpotTotal{i,j}=sum(IcfpM{i}{j},2);
        IrfpSpotTotal{i,j}=sum(IrfpM{i}{j},2);
        
        %find number of spots
        %find nonzero elements
        [RowYfpnz,ColYfpnz]=find(IyfpM{i}{j}>0);
        %find multitude of spots by sorting and counting occurences by
        %sorting
        [MultiSpotYFP,SortIdxYfp]=sort(RowYfpnz);         
        %sort, occurrence, then unique and link
        DD=MultiSpotYFP;
        [aanz,bbnz]=hist(DD,unique(DD)); % aa is multitude, bb is row number
        
        % now we now each row's multitude, thus per frame the number of spots. 
        % Per cell j we need to create one array with spot numbers from
        % first to second frame. So we need to fill in the gaps because we
        % only have the nonzero elements now 
        YFPSpotNumber{i,j}=zeros(1,frames{i}{j});
        for L=1:size(bbnz,1)
        YFPSpotNumber{i,j}(bbnz(L))=aanz(L);
        end
        
        % do same for RFP
        [RowRfpnz,ColRfpnz]=find(IrfpM{i}{j}>0);
        [MultiSpotRFP,SortIdxRfp]=sort(RowRfpnz);         
        DD=MultiSpotRFP;
        [aanz,bbnz]=hist(DD,unique(DD)); 
        
        RFPSpotNumber{i,j}=zeros(1,frames{i}{j});

        for L=1:size(bbnz,1)
        RFPSpotNumber{i,j}(bbnz(L))=aanz(L);
        end
    end
end 



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
colormap(fig4,jet) % heat map
xlabel('Cell Length'); ylabel('Position in Cell');
title('Kymo data: CFP');
grid on

figure(5)
subplot(1,3,1)
hold on
hist3(Numcfp','Edges',thisedge2)
colormap(fig5,jet) % heat map
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

%Filter on Intensity
IYFP=Iyfp*Intensityval(2);
[bin_yfp,idx_yfp]=find(IYFP>1100); %980=Noise-Fraction 0.57

NumYFP(:,idx_yfp)=Numyfp(:,idx_yfp);
NumYFPnz(1,:)=nonzeros(NumYFP(1,:));
NumYFPnz(2,:)=nonzeros(NumYFP(2,:));

FilteredSignalYFP=size(NumYFPnz,2)/size(IYFP,2); % Percentage of total signal

figure(4)
subplot(1,3,2)
Heatmap = hist3(NumYFPnz','Edges',thisedge2);
pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
xlabel('Cell Length'); ylabel('Position in Cell');
title('Kymo data: YFP');
grid on

figure(5)
subplot(1,3,2)
hold on
hist3(Numyfp','Edges',thisedge2)
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
xlabel('Cell Length'); ylabel('Position in Cell');
title('Kymo data: RFP');
grid on

figure(5)
subplot(1,3,3)
hold on
hist3(Numrfp','Edges',thisedge2)
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
YFPcal=953;
RFPcal=990;

% CFP

plotcfp(1,:) = Lcfp;
plotcfp(2,:) = Fcfp;
plotcfp(3,:) = Icfp*Intensityval(1);
plotcfp = unique(plotcfp','rows')';

subplot(1,3,1)
hold on
scatter(plotcfp(1,:),plotcfp(2,:),'b','x');
scatter(plotcfp(1,:),plotcfp(3,:),'r','x');
myfit=polyfit(plotcfp(1,:),plotcfp(2,:),1);
x=0:0.01:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Cell Time'); ylabel('Normalized full cell intensity'); 
title('Kymo data: CFP')
hold off
% axis([0 1 -0.1 1.1])

% YFP

plotyfp(1,:) = Lyfp;
plotyfp(2,:) = Fyfp/YFPcal;
plotyfp(3,:) = Iyfp*Intensityval(2)/YFPcal;
plotyfp = unique(plotyfp','rows')';

subplot(1,3,2)
hold on
scatter(plotyfp(1,:),plotyfp(2,:),'b','x');
% scatter(plotyfp(1,:),plotyfp(3,:),'r','x');
for i=channels
    for j=1:Ncells{i}
scatter((1:size(IyfpSpotTotal{i,j},1))/size(IyfpSpotTotal{i,j},1),IyfpSpotTotal{i,j}./YFP.filterval,'r')
    end
end
myfit=polyfit(plotyfp(1,:),plotyfp(2,:),1);
x=0:0.01:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Normalized Cell Time (-)'); ylabel('Number of Tus (-)'); 
title('Tus Stoichiometry vs. Time')
hold off
axis([0 1 0.01 30])
set(gca,'FontSize',18)


% RFP

plotrfp(1,:) = Lrfp;
plotrfp(2,:) = Frfp;
plotrfp(3,:) = Irfp*Intensityval(3);
plotrfp = unique(plotrfp','rows')';

subplot(1,3,3)
hold on
scatter(plotrfp(1,:),plotrfp(2,:),'b','x');
% scatter(plotrfp(1,:),plotrfp(3,:),'r','x');
for i=channels
    for j=1:Ncells{i}
scatter((1:size(IrfpSpotTotal{i,j},1))/size(IrfpSpotTotal{i,j},1),IrfpSpotTotal{i,j},'r')
    end
end
myfit=polyfit(plotrfp(1,:),plotrfp(2,:),1);
x=0:0.01:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Cell Time'); ylabel('Normalized full cell intensity'); 
title('Kymo data: RFP')
hold off
axis([0 1 10 70000])

%% YFP Spot number vs time
fig7=figure(7);

YfpSpot=[];
Frame=[];

for i=channels
    for j=1:Ncells{i}
        YfpSpot=[YfpSpot YFPSpotNumber{i,j}];
        Frame=[Frame (1:frames{i}{j})/frames{i}{j}];
    end
end

for i=1:max(YfpSpot);
p{i} = find(YfpSpot == i);
m(i) = mean(Frame(p{i}));
Spotstd(i)=std(Frame(p{i}));
end

Y=linspace(0,max(YfpSpot),max(YfpSpot));

Xl=m-Spotstd;
Xr=m+Spotstd;

hold on
scatter(Frame,YfpSpot,'r','filled')
plot(m,Y,'b--','LineWidth',3)
plot(Xl,Y,'b','LineWidth',1)
plot(Xr,Y,'b','LineWidth',1)
hold off
title('Tus Spot number vs. Normalized Cell cycle')
xlabel('Normalized Cell Time (-)'); ylabel('Number of Spots (-)')
set(gca,'FontSize',18)




%%
hold on

for i=channels
    for j=1:Ncells{i}
        
scatter((1:size(IrfpSpotTotal{i,j},1))/size(IrfpSpotTotal{i,j},1),IrfpSpotTotal{i,j}/RFP.filterval,'r','filled')
         
    end
end
