clear all
close
clc

folder='C:\Users\water\Documents\GitHub\Data\Target Data\Kymo Data\Results';
slash = '\';
channels = 2;
Intensityval = [500,250,100];

Acfp=[];    Ayfp=[];    Arfp=[];
Bcfp=[];    Byfp=[];    Brfp=[];
Ccfp=[];    Cyfp=[];    Crfp=[];

for i=1:channels;
    E{i}=load(strcat(folder,slash,'Results_Ch',num2str(i),'.mat')); 
end

allCFP_L = [];


for i=1:channels
    
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
%% CFP

figure(1)
hold on
scatter(single(Acfp),Bcfp,Ccfp,'b','filled');
% myfit=polyfit(Acfp,Bcfp,4);
% x=0:0.01:1;
% y=polyval(myfit,x);
% plot(x,y,'r','LineWidth',5)
hold off
axis([0 1 0 1])


%% YFP

figure(2)
hold on
scatter(single(Ayfp),Byfp,Cyfp,'m','filled');
% myfit=polyfit(Ayfp,Byfp,4);
% x=0:0.01:1;
% y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
hold off
axis([0 1 0 1])

%% RFP

figure(3)
hold on
scatter(single(Arfp),Brfp,Crfp,'r','filled');
% myfit=polyfit(Arfp,Brfp,4);
% x=0:0.01:1;
% y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
hold off
axis([0 1 0 1])
        
        