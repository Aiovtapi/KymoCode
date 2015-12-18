% Code to normalise intensity data to compare different data sets for the
% doubling event of the intensity. This has been created since initial
% intensities can differ significantly, so that we cannot compare two sets
% of data in terms of doubling. Run GaussCalcs before.
clear all
close all

GaussCalcs

% 4 is a single bac instead of two 

% Kdi(:,4)=2*Kdi(:,4); KdiR1(:,4)=2*KdiR1(:,4); 
% KdiR2(:,4)=2*KdiR2(:,4); KdiR3(:,4)=2*KdiR3(:,4);

KdiTotal=Kdi+KdiR1+KdiR2+KdiR3;
KdiNorm=zeros(MaxBacLifed,Ncells);
KdiRatio=KdiTotal./KdiFC;

for i=1:MaxBacLifed
    for j=1:Ncells
        KdiNorm(i,j)=KdiTotal(i,j)/(mean(KdiTotal(1:5,j)));
    end
end
  
for i=1:MaxBacLifed
    MdiNorm(i,1)=mean(nonzeros(KdiNorm(i,:)));
    
    MdiNorm(i,2)=std(nonzeros(KdiNorm(i,:)));
    
    MdiNormFC(i,1)=mean(KdiRatio(i,:));
    MdiNormFC(i,2)=std(KdiRatio(i,:));
end


X=linspace(0,1,length(KdiNorm(:,1)));

%% Plots
figure(2)
hold on 
plot(X,Md(:,1),'-ob','LineWidth',3);
hold off
title('R2 Spot Intensity vs. Time','FontSize',24);
xlabel('Normalised Cell Time (-)','FontSize',18,'FontWeight','bold');
ylabel('Normalized Intensity (-)','FontSize',18,'FontWeight','bold');
txtbox1={'Ncells=40'};
annotation('textbox', [0.15,0.8,0.1,0.1],...
           'String', txtbox1,'FontSize',20,'FontWeight','bold')
set(gca,'FontSize',20);
%%
figure(3)
hold on 
plot(X,(MdiNorm(:,1)),'-or','LineWidth',3);
hold off
title('Normalised R2 Spot Intensity vs. Time','FontSize',24);
xlabel('Normalised Cell Time (-)','FontSize',18,'FontWeight','bold');
ylabel('Normalized Intensity (-)','FontSize',18,'FontWeight','bold');
txtbox1={'Ncells=40'};
line([0 1],[1.8 1.8]); line([0 1],[2 2]);
axis([0 1 1 3])
annotation('textbox', [0.15,0.8,0.1,0.1],...
           'String', txtbox1,'FontSize',20,'FontWeight','bold')
set(gca,'FontSize',20);
%% Plots
figure(4)
hold on 
plot(X,M(:,1)*2/3600,'-ob','LineWidth',3);
hold off
title('Replisome DnaN Proteins vs. Time','FontSize',24);
xlabel('Normalised Cell Time (-)','FontSize',18,'FontWeight','bold');
ylabel('Number (-)','FontSize',18,'FontWeight','bold');
txtbox1={'Ncells=8'};
annotation('textbox', [0.15,0.8,0.1,0.1],...
           'String', txtbox1,'FontSize',20,'FontWeight','bold')
set(gca,'FontSize',20);

%% 
figure(5)
hold on
plot(X,MdiDifX.^2,'-ok','LineWidth',3)
hold off
title('Normalised Distance Between Spots vs. Time','FontSize',24);
xlabel('Normalised Cell Time (-)','FontSize',18,'FontWeight','bold');
ylabel('Number (-)','FontSize',18,'FontWeight','bold');
txtbox1={'Ncells=40'};
annotation('textbox', [0.15,0.8,0.1,0.1],...
           'String', txtbox1,'FontSize',20,'FontWeight','bold')
set(gca,'FontSize',20);