clear all
close
clc

folder='/Users/rleeuw/Work/Data/160205_BN2384_and_Beam_Profiles';
exps=[1 4 5 7];

Acfp=[];    Ayfp=[];    Arfp=[];
Bcfp=[];    Byfp=[];    Brfp=[];
Ccfp=[];    Cyfp=[];    Crfp=[];

    j=1;
    
for i=exps;
    E{j}=load(strcat(folder,'/',num2str(i),'/Results.mat')); 
    imflip{j}=load(strcat(folder,'/',num2str(i),'/imgflip.mat'));
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
    CFP{i,j}{1}=[];
else
        for k=1:NspotsCFP
                if iscell(CFPld{i,j}{k})
                    
                    CFP{i,j}{k}=CFPld{i,j}{k}{k};
                    if LNormCFP{i,j}(1)>0.5
                        CFP{i,j}{k}=CellLength{i,j}-CFP{i,j}{k};
                    end
                else
                    
                    CFP{i,j}{k}=CFPld{i,j}{k}; 
                    if LNormCFP{i,j}(1)>0.5
                        CFP{i,j}{k}=CellLength{i,j}-CFP{i,j}{k};
                    end
                end
                
                if ~iscell(CFPld{i,j}{k})
                Acfp=[Acfp CellLength{i,j}];
                Bcfp=[Bcfp CFP{i,j}{k}(1,2)/CellLength{i,j}];
                Ccfp=[Ccfp CFPld{i,j}{k}(1,1)/700];
                end
        end
            
end
        
if NspotsYFP==0
    YFP{i,j}{1}=[];
else
        for k=1:NspotsYFP
                if iscell(YFPld{i,j}{k})
                    YFP{i,j}{k}=YFPld{i,j}{k}{k};
                    if LNormCFP{i,j}(1)>0.5 && ~iscell(YFP{i,j}{k})
                        YFP{i,j}{k}=CellLength{i,j}-YFP{i,j}{k};
                    end
                else
                    YFP{i,j}{k}=YFPld{i,j}{k}; 
                        if LNormCFP{i,j}(1)>0.5
                            YFP{i,j}{k}=CellLength{i,j}-YFP{i,j}{k};
                        end
                end
                
                if ~iscell(YFPld{i,j}{k})
                Ayfp=[Ayfp CellLength{i,j}];
                Byfp=[Byfp YFP{i,j}{k}(1,2)/CellLength{i,j}];
                Cyfp=[Cyfp YFPld{i,j}{k}(1,1)/300];
                end
        end
            
end

if NspotsRFP==0
    RFP{i,j}{1}=[];
else
        for k=1:NspotsRFP
                if iscell(RFPld{i,j}{k})
                    
                    RFP{i,j}{k}=RFPld{i,j}{k}{k};
                    if LNormCFP{i,j}(1)>0.5
                        RFP{i,j}{k}=CellLength{i,j}-RFP{i,j}{k};
                    end
                else
                    
                    RFP{i,j}{k}=RFPld{i,j}{k}; 
                    if LNormCFP{i,j}(1)>0.5
                        RFP{i,j}{k}=CellLength{i,j}-RFP{i,j}{k};
                    end
                end
                
                if ~iscell(RFPld{i,j}{k})
                Arfp=[Arfp CellLength{i,j}];
                Brfp=[Brfp RFP{i,j}{k}(1,2)/CellLength{i,j}];
                Crfp=[Crfp RFPld{i,j}{k}(1,1)/700];
                end
        end
            
end

%         for k=1:NspotsYFP
%             if iscell(YFPld{i,j}{k});
%                YFP{i,j}{k}=YFPld{i,j}{k}{k};
%             else
%                YFP{i,j}{k}=YFPld{i,j}{k}; 
%             end
%         end
%         
%         for k=1:NspotsRFP
%             if iscell(RFPld{i,j}{k});
%                RFP{i,j}{k}=RFPld{i,j}{k}{k};
%             else
%                RFP{i,j}{k}=RFPld{i,j}{k}; 
%             end
%         end

    end
end
    


%% CFP

figure(1)
hold on
scatter(single(Acfp),Bcfp,Ccfp,'b','filled');
myfit=polyfit(Acfp,Bcfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
plot(x,y,'r','LineWidth',5)
hold off
axis([12 43 -0.1 1.1])


%% YFP

figure(2)
hold on
scatter(single(Ayfp),Byfp,Cyfp,'m','filled');
myfit=polyfit(Ayfp,Byfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
hold off
axis([12 43 -0.1 1.1])
axis([12 43 -0.1 1.1])

%% RFP

figure(3)
hold on
scatter(single(Arfp),Brfp,Crfp,'r','filled');
myfit=polyfit(Arfp,Brfp,4);
x=15:0.1:45;
y=polyval(myfit,x);
% plot(x,y,'k','LineWidth',5)
hold off
axis([12 43 -0.1 1.1])
axis([12 43 -0.1 1.1])

        
        