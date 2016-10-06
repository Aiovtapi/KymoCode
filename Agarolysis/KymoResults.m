clear 
close all
clc

folder='C:\Users\water\Documents\GitHub\Data\Target Data\Kymo Data\Results';
slash = '\';
channels = [1,2];
Intensityval = [350,60,200];

Acfp=[];    Ayfp=[];    Arfp=[];
Bcfp=[];    Byfp=[];    Brfp=[];
Ccfp=[];    Cyfp=[];    Crfp=[];

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
                    Acfp=[Acfp h/frames{i}{j}];
                    Bcfp=[Bcfp CFPx{i}{j}{k}(h,2)/length];
                    Ccfp=[Ccfp CFPx{i}{j}{k}(h,1)/Intensityval(1)];
                end
            end
        end
        
        if NspotsYFP==0
            YFPx{i}{j}{1}=[];
        else
            for k=1:NspotsYFP
                for h = 1:frames{i}{j};
                    length = size(E{i}.DataStruct(2,j).ydatacrpdR1{h,k},2);
                    Ayfp=[Ayfp h/frames{i}{j}];
                    Byfp=[Byfp YFPx{i}{j}{k}(h,2)/length];
                    Cyfp=[Cyfp YFPx{i}{j}{k}(h,1)/Intensityval(2)];
                end
            end
        end
        
        if NspotsRFP==0
            RFPx{i}{j}{1}=[];
        else
            for k=1:NspotsRFP
                for h = 1:frames{i}{j};
                    length = size(E{i}.DataStruct(3,j).ydatacrpdR1{h,k},2);
                    Arfp=[Arfp h/frames{i}{j}];
                    Brfp=[Brfp RFPx{i}{j}{k}(h,2)/length];
                    Crfp=[Crfp RFPx{i}{j}{k}(h,1)/Intensityval(3)];
                end
            end
        end
    end
end

Yremove = Cyfp == 0;
Cyfp(Yremove) = [];
Ayfp(Yremove) = [];
Byfp(Yremove) = [];

Rremove = Crfp == 0;
Crfp(Rremove) = [];
Arfp(Rremove) = [];
Brfp(Rremove) = [];

Cremove = Ccfp == 0;
Ccfp(Cremove) = [];
Acfp(Cremove) = [];
Bcfp(Cremove) = [];

%% Position vs. cell length
% CFP

fig1 = figure(1);
set(fig1,'Position',[20,300,1800,500])
subplot(1,3,1)
hold on
scatter(single(Acfp),Bcfp,Ccfp,'b','filled');
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
scatter(single(Ayfp),Byfp,Cyfp,'m','filled');
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
scatter(single(Arfp),Brfp,Crfp,'r','filled');
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

fig1 = figure(2);
set(fig1,'Position',[20,300,1800,500])
subplot(1,3,1)
hold on
scatter(Bcfp,Ccfp,'b','x');
myfit=polyfit(Bcfp,Ccfp,4);
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
scatter(Byfp,Cyfp,'m','x');
myfit=polyfit(Byfp,Cyfp,4);
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
scatter(Brfp,Crfp,'r','x');
myfit=polyfit(Brfp,Crfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Kymo data: RFP')
hold off
axis([0 1 -0.1 35])

%% Numspots vs. position
% CFP
fig1 = figure(3);
set(fig1,'Position',[20,300,1800,500])
subplot(1,3,1)

[numbin,edges] = histcounts(Bcfp,20);
norm = max(numbin)/35;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Bcfp,Ccfp,'b','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Kymo data: CFP')
axis([0 1 -0.1 40])

% YFP
figure(3);
subplot(1,3,2)

[numbin,edges] = histcounts(Byfp,20);
norm = max(numbin)/35;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Byfp,Cyfp,'m','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Kymo data: YFP')
axis([0 1 -0.1 40])

% RFP
figure(3);
subplot(1,3,3)

[numbin,edges] = histcounts(Brfp,20);
norm = max(numbin)/35;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Brfp,Crfp,'r','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)'); 
title('Kymo data: RFP')
axis([0 1 -0.1 40])