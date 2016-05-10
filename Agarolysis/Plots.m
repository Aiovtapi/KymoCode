hold on
scatter(CellLength(find(singletus)),nonzeros(singletus),'r','filled')
scatter(CellLength(find(doubletus(:,1))),nonzeros(doubletus(:,1)),'b','filled')
scatter(CellLength(find(doubletus(:,2))),nonzeros(doubletus(:,2)),'b','filled')
scatter(CellLength(find(tripletus(:,1))),nonzeros(tripletus(:,1)),'k','filled')
scatter(CellLength(find(tripletus(:,2))),nonzeros(tripletus(:,2)),'k','filled')
scatter(CellLength(find(tripletus(:,3))),nonzeros(tripletus(:,3)),'k','filled')
scatter(CellLength(find(quadtus(:,1))),nonzeros(quadtus(:,1)),'m','filled')
scatter(CellLength(find(quadtus(:,1))),nonzeros(quadtus(:,2)),'m','filled')
scatter(CellLength(find(quadtus(:,1))),nonzeros(quadtus(:,3)),'m','filled')
scatter(CellLength(find(quadtus(:,1))),nonzeros(quadtus(:,4)),'m','filled')
scatter(CellLength(find(quinttus(:,1))),nonzeros(quinttus(:,1)),'c','filled')
scatter(CellLength(find(quinttus(:,1))),nonzeros(quinttus(:,2)),'c','filled')
scatter(CellLength(find(quinttus(:,1))),nonzeros(quinttus(:,3)),'c','filled')
scatter(CellLength(find(quinttus(:,1))),nonzeros(quinttus(:,4)),'c','filled')
scatter(CellLength(find(quinttus(:,1))),nonzeros(quinttus(:,5)),'c','filled')
hold off

%% Tus spot Numbers histogram by length
figure(1)
axis([0 6 0.01 0.7])
subplot(1,3,1);
h1=histogram(nonzeros(TusSpotNumber(find(CellLength>10 & CellLength<20))),'FaceColor','r','facealpha',.5,'Normalization','probability');
set(gca,'XTick',[0 1 2 3 4 5 6])
set(gca,'FontSize',12)
axis([0 6 0 0.5])
xlabel('Spots (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Tus spot number - 1.6 to 3.2 ?m')
subplot(1,3,2);
h2=histogram(nonzeros(TusSpotNumber(find(CellLength>20 & CellLength<30))),'FaceColor','b','facealpha',.5,'Normalization','probability');
set(gca,'XTick',[0 1 2 3 4 5 6])
set(gca,'FontSize',12)
axis([0 6 0 0.5])
xlabel('Spots (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Tus spot number - 3.2 to 4.8 ?m')
subplot(1,3,3);
h3=histogram(nonzeros(TusSpotNumber(find(CellLength>30 & CellLength<40))),'FaceColor','y','facealpha',.5,'Normalization','probability');
set(gca,'XTick',[0 1 2 3 4 5 6])
set(gca,'FontSize',12)
axis([0 6 0 0.5])
xlabel('Spots (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Tus spot number - 4.8 to 6.4 ?m')

%% dif spot numbers
figure(2)
hold on
axis([0 4 0.01 1.1])
h1=histogram(nonzeros(difSpotNumber(find(CellLength>10 & CellLength<20))),'BinWidth',.9,'FaceColor','r','facealpha',.5,'Normalization','probability');
h2=histogram(nonzeros(difSpotNumber(find(CellLength>20 & CellLength<30))),'BinWidth',.9,'FaceColor','b','facealpha',.5,'Normalization','probability');
h3=histogram(nonzeros(difSpotNumber(find(CellLength>30 & CellLength<40))),'BinWidth',.9,'FaceColor','y','facealpha',.5,'Normalization','probability');
h4=histogram(nonzeros(difSpotNumber(find(CellLength>40 & CellLength<50))),'BinWidth',.9,'FaceColor','m','facealpha',.5,'Normalization','probability');
hold off
set(gca,'XTick',[0 1 2 3 4])
set(gca,'FontSize',16)
xlabel('Spots (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('R2 spots vs. cell length')
legend('1.6 to 3.2 ?m','3.2 to 4.8 ?m','4.8 to 6.4 ?m')

%% Distance plot
smallT=NormTusPos(find(CellLength>10 & CellLength<20));
mediumT=NormTusPos(find(CellLength>20 & CellLength<30));
largeT=NormTusPos(find(CellLength>30 & CellLength<40));
allT=NormTusPos(find(CellLength>10 & CellLength<40));

CLAllVec=CellLength(find(CellLength>10 & CellLength<40));

smalld=NormdifPos(find(CellLength>10 & CellLength<20));
mediumd=NormdifPos(find(CellLength>20 & CellLength<30));
larged=NormdifPos(find(CellLength>30 & CellLength<40));
alld=NormdifPos(find(CellLength>10 & CellLength<40));

for i=1:length(allT)
    
        % 1 Tus
    if (size(allT{i},2)==1 && size(alld{i},2)>0);
        sd{i}=min(abs(alld{i}-allT{i}));
        
        % 2 Tus
    elseif size(allT{i},2)==2 && size(alld{i},2)>0;
        sd{i}(1)=min(abs(alld{i}-allT{i}(1)));
        sd{i}(2)=min(abs(alld{i}-allT{i}(2)));
        
        % 3 Tus
    elseif size(allT{i},2)==3 && size(alld{i},2)>0;
        sd{i}(1)=min(abs(alld{i}-allT{i}(1)));
        sd{i}(2)=min(abs(alld{i}-allT{i}(2)));
        sd{i}(3)=min(abs(alld{i}-allT{i}(3)));
        
        % 4 Tus
    elseif size(allT{i},2)==4 && size(alld{i},2)>0;
        sd{i}(1)=min(abs(alld{i}-allT{i}(1)));
        sd{i}(2)=min(abs(alld{i}-allT{i}(2)));
        sd{i}(3)=min(abs(alld{i}-allT{i}(3)));
        sd{i}(4)=min(abs(alld{i}-allT{i}(4)));
        
        % 5 Tus
    elseif size(allT{i},2)==5 && size(alld{i},2)>0;
        sd{i}(1)=min(abs(alld{i}-allT{i}(1)));
        sd{i}(2)=min(abs(alld{i}-allT{i}(2)));
        sd{i}(3)=min(abs(alld{i}-allT{i}(3)));
        sd{i}(4)=min(abs(alld{i}-allT{i}(4)));
        sd{i}(5)=min(abs(alld{i}-allT{i}(5)));
        
    else
        sd{i}=0;
        
    end
    sizesd(i)=size(sd{i},2);
end

%%
for i=1:length(sizesd)
    if sizesd(i)==1
    singlespotdist(i)=sd{i};
    elseif sizesd(i)==2
        doublespotdist(i,:)=sd{i};
    elseif sizesd(i)==3
        triplespotdist(i,:)=sd{i};
    elseif sizesd(i)==4
        quadspotdist(i,:)=sd{i};
    elseif sizesd(i)==5
        quintspotdist(i,:)=sd{i};
    else
        singlespotdist(i)=0;
        doublespotdist(i,1)=0; doublespotdist(i,2)=0;
        triplespotdist(i,1:3)=0;
        continue
    end
end
hold on
 scatter(CLAllVec,[singlespotdist 0],'b','filled')
%  scatter(CLAllVec,doublespotdist(:,1),'r','filled')
%  scatter(CLAllVec,doublespotdist(:,2),'r','filled')
% scatter(CLAllVec,triplespotdist(:,1),'k','filled')
% scatter(CLAllVec,triplespotdist(:,2),'k','filled')
% scatter(CLAllVec,triplespotdist(:,3),'k','filled')
axis([15 40 0.01 0.6])

%% Hist Dist
L3=length(nonzeros(singlespotdist()));
figure(3)
hold on
h1=histogram(nonzeros(singlespotdist),'BinWidth',.05,'FaceColor','r','facealpha',.5,'Normalization','probability');
% histfit(nonzeros(singlespotdist),10,'exponential')
axis([0 .6 0 0.5])
set(gca,'XTick',[0 .1 .2 .3 .4 .5 .6 .7 .8 1])
set(gca,'FontSize',16)
xlabel('Normalised Distance (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Distance between single Tus spot and R2 spot N=26');

%%
figure(4)
L4=length(nonzeros(doublespotdist(find(CLAllVec<50))));
h2=histogram(nonzeros(doublespotdist(find(CLAllVec<50))),10,'BinWidth',.05,'FaceColor','b','facealpha',.5,'Normalization','probability');
axis([0 0.7 0 .2])
set(gca,'XTick',[0 .1 .2 .3 .4 .5 .6 .7 .8 1])
set(gca,'FontSize',16)
xlabel('Normalised Distance (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Distance between double Tus spots and R2 spot N=23');
%%
figure(5)
L5=length(nonzeros(triplespotdist(find(CLAllVec<50))));
h3=histogram(nonzeros(triplespotdist(find(CLAllVec<50))),10,'BinWidth',.07,'FaceColor','g','facealpha',.5,'Normalization','probability');
axis([0 .7 0 .35])
set(gca,'XTick',[0 .1 .2 .3 .4 .5 .6 .7 .8 1])
set(gca,'FontSize',16)
xlabel('Normalised Distance (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Distance between triple Tus spot and R2 spot N=16');
%%
figure(6)
L6=length(nonzeros(quadspotdist(find(CLAllVec<50))));
h4=histogram(nonzeros(quadspotdist(find(CLAllVec<50))),'BinWidth',.05,'FaceColor','g','facealpha',.5,'Normalization','probability');
axis([0 .6 0 .5])
set(gca,'XTick',[0 .1 .2 .3 .4 .5 .6 .7 .8 1])
set(gca,'FontSize',16)
xlabel('Normalised Distance (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Distance between quad Tus spot and R2 spot N=7');
%%
figure(7)
L7=length(nonzeros(quintspotdist(find(CLAllVec<50))));
h5=histogram(nonzeros(quintspotdist(find(CLAllVec<50))),10,'BinWidth',.06,'FaceColor','g','facealpha',.5,'Normalization','probability');
axis([0 .7 0 .35])
set(gca,'XTick',[0 .1 .2 .3 .4 .5 .6 .7 .8 1])
set(gca,'FontSize',16)
xlabel('Normalised Distance (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Distance between five Tus spots and R2 spot N=16');