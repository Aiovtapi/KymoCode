clear all
close
clc

folder='C:\Users\water\Documents\GitHub\Data\Target Data\Agar Data';
slash = '\';
exps=[1 2  3 4 5 7 8 9];
Intensityval = [700, 300, 700];

Acfp=[];    Ayfp=[];    Arfp=[];
Bcfp=[];    Byfp=[];    Brfp=[];
Ccfp=[];    Cyfp=[];    Crfp=[];


    
j=1; 
for i=exps;
    E{j}=load(strcat(folder,slash,num2str(i),slash,'Results.mat')); 
%     imflip{j}=load(strcat(folder,slash,num2str(i),slash,'imgflip.mat'));
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
                Acfp=[Acfp CellLength{i,j}];
                Bcfp=[Bcfp CFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Ccfp=[Ccfp CFPld{i,j}{k}(1,1)/Intensityval(1)];
            end
        end
        
        if NspotsYFP==0
            YFPld{i,j}{1}=[];
        else
            for k=1:NspotsYFP
                Ayfp=[Ayfp CellLength{i,j}];
                Byfp=[Byfp YFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Cyfp=[Cyfp YFPld{i,j}{k}(1,1)/Intensityval(2)];
            end
        end
        
        if NspotsRFP==0
            RFPld{i,j}{1}=[];
        else
            for k=1:NspotsRFP            
                Arfp=[Arfp CellLength{i,j}];
                Brfp=[Brfp RFPld{i,j}{k}(1,2)/CellLength{i,j}];
                Crfp=[Crfp RFPld{i,j}{k}(1,1)/Intensityval(3)];
            end
        end
    end
end

remove = Crfp<70;
Cyfp(remove) = [];
Ayfp(remove) = [];
Byfp(remove) = [];

%% Position vs. cell length
% CFP

fig1 = figure(1);
set(fig1,'Position',[20,300,1800,500])
subplot(1,3,1)

figure(1)
hold on
scatter(single(Acfp),Bcfp,Ccfp,'b','filled');
myfit=polyfit(Acfp,Bcfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
plot(x,y,'r','LineWidth',5)
xlabel('Cell Length'); ylabel('Normalized osition of spot in cell'); 
title('Agar data: CFP')
hold off
axis([12 43 -0.1 1.1])


% YFP

subplot(1,3,2);
hold on
scatter(single(Ayfp),Byfp,Cyfp,'m','filled');
myfit=polyfit(Ayfp,Byfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
xlabel('Cell Length'); ylabel('Normalized osition of spot in cell'); 
title('Agar data: YFP')
hold off
axis([12 43 -0.1 1.1])
axis([12 43 -0.1 1.1])


% RFP

subplot(1,3,3)
hold on
scatter(single(Arfp),Brfp,Crfp,'r','filled');
myfit=polyfit(Arfp,Brfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
xlabel('Cell Length'); ylabel('Normalized osition of spot in cell'); 
title('Agar data: RFP')
% plot(x,y,'k','LineWidth',5)
hold off
axis([12 43 -0.1 1.1])
axis([12 43 -0.1 1.1])


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
title('Agar data: CFP')
hold off
axis([0 1 -0.1 90])


% YFP

subplot(1,3,2)
hold on
scatter(Byfp,Cyfp,'m','x');
myfit=polyfit(Byfp,Cyfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Agar data: YFP')
hold off
axis([0 1 -0.1 90])


% RFP

subplot(1,3,3)
hold on
scatter(Brfp,Crfp,'r','x');
myfit=polyfit(Brfp,Crfp,4);
x=0:00.1:1;
y=polyval(myfit,x);
plot(x,y,'k','LineWidth',3)
xlabel('Position in cell'); ylabel('Spot Intensity'); 
title('Agar data: RFP')
hold off
axis([0 1 -0.1 90])


%% Numspots vs. position
% CFP
bins = 15;
thisedge = (0:bins)/bins;

fig1 = figure(3);
set(fig1,'Position',[20,300,1800,500])
subplot(1,3,1)

[numbin,edges] = histcounts(Bcfp,thisedge);
norm = max(numbin)/80;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Bcfp,Ccfp,'b','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Agar data: CFP')
axis([0 1 -0.1 90])

% YFP
figure(3);
subplot(1,3,2)

[numbin,edges] = histcounts(Byfp,thisedge);
norm = max(numbin)/80;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Byfp,Cyfp,'m','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)');
title('Agar data: YFP')
axis([0 1 -0.1 90])

% RFP
figure(3);
subplot(1,3,3)

[numbin,edges] = histcounts(Brfp,thisedge);
norm = max(numbin)/80;
X = diff(edges);
X = cumsum(X) - X(1)/2;
hold on
scatter(Brfp,Crfp,'r','x');
plot(X,numbin/norm,'k','LineWidth',3)
hold off
xlabel('Position in cell'); ylabel('Amount of spots (normalized)'); 
title('Agar data: RFP')
axis([0 1 -0.1 90])
        