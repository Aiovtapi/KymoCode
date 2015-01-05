load Foci_WT_DnaN_Avg.txt -ascii;   
load Foci_WT_DnaN_HIGH.txt -ascii;
load Foci_WT_DnaN_LOW.txt -ascii;


load Foci_GammaMinus_DnaN_Avg.txt -ascii;   
load Foci_GammaMinus_DnaN_HIGH.txt -ascii;
load Foci_GammaMinus_DnaN_LOW.txt -ascii;

load Hist_Data_GammaMinus.txt -ascii;
load Hist_Data_WT.txt -ascii;


load CytoplasmicConcentration_WT_AVG.txt -ascii;
load CytoplasmicConcentration_WT_LOW.txt -ascii;
load CytoplasmicConcentration_WT_HIGH.txt -ascii;

load CytoplasmicConcentration_GammaMinus_AVG.txt -ascii;
load CytoplasmicConcentration_GammaMinus_LOW.txt -ascii;
load CytoplasmicConcentration_GammaMinus_HIGH.txt -ascii;



load Cyto_WT_DnaN_Avg.txt -ascii;
load Cyto_WT_DnaN_LOW.txt -ascii;  
load Cyto_WT_DnaN_HIGH.txt -ascii;

WT_X_Cyto = CytoplasmicConcentration_WT_AVG(:,1);
WT_AVG_Cyto = CytoplasmicConcentration_WT_AVG(:,2);
WT_hi_Cyto = CytoplasmicConcentration_WT_HIGH(:,2);
WT_lo_Cyto = CytoplasmicConcentration_WT_LOW(:,2);

GammaMinus_X_Cyto = CytoplasmicConcentration_GammaMinus_AVG(:,1);
GammaMinus_AVG_Cyto = CytoplasmicConcentration_GammaMinus_AVG(:,2);
GammaMinus_hi_Cyto = CytoplasmicConcentration_GammaMinus_HIGH(:,2);
GammaMinus_lo_Cyto = CytoplasmicConcentration_GammaMinus_LOW(:,2);

gammaNegX_hist = Hist_Data_GammaMinus(:,1);
gammaNegY_hist =Hist_Data_GammaMinus(:,2);

WT_X_hist = Hist_Data_WT(:,1);
WT_Y_hist = Hist_Data_WT(:,2);

X_WT = Foci_WT_DnaN_Avg(:,1);     
WT_AVG = Foci_WT_DnaN_Avg(:,2);  
WT_lo = Foci_WT_DnaN_LOW(:,2);
WT_hi = Foci_WT_DnaN_HIGH(:,2);


X_GammaMinus = Foci_GammaMinus_DnaN_Avg(:,1);     
GammaMinus_AVG = Foci_GammaMinus_DnaN_Avg(:,2);  
GammaMinus_lo = Foci_GammaMinus_DnaN_LOW(:,2);
GammaMinus_hi = Foci_GammaMinus_DnaN_HIGH(:,2);



%%Plotting the foci content
h=figure;
ciplot(WT_lo,WT_hi,X_WT,[51/255 153/255 1]);
hold on;
%ciplot(GammaMinus_lo,GammaMinus_hi,X_GammaMinus,[1 153/255 153/255]);   

% 
plot(X_WT,WT_AVG,'-b','LineWidth',4); hold on;

%plot(X_GammaMinus,GammaMinus_AVG,'-r','LineWidth',4); hold on;  

 set(gca,'TickLength',[0.02 0.02]);
   set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
   box off;
   axis([-0.25 1.25 0 70])
   
   ylabel('\beta_2 clamps (-)', 'fontsize', 26, 'fontweight', 'bold')
   xlabel('Normalized time (-)', 'fontsize', 26, 'fontweight', 'bold')
   
      
  %  hleg = legend('YPet-DnaN \gamma+','YPet-DnaN \gamma-');

   set(hleg, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
   
     set(gca, 'YTick', [0 10 20 30 40 50 60 70]);
    set(gca, 'XTick', [-0.25 0 0.25 0.5 0.75 1 1.25]);
%            

%print(h, '-dpdf', '-r600','GammaPos_GammaNeg_CoPlot')


hold off;


















% h3=figure;
% %WT_X_Cyto = CytoplasmicConcentration_WT_AVG(:,1);
% %WT_AVG_Cyto = CytoplasmicConcentration_WT_AVG(:,2);
% ciplot(WT_lo_Cyto,WT_hi_Cyto,X_WT,[51/255 153/255 1]);
% hold on;
% ciplot(GammaMinus_lo_Cyto,GammaMinus_hi_Cyto,GammaMinus_X_Cyto,[1 153/255 153/255]);  
% plot(WT_X_Cyto,WT_AVG_Cyto,'-b','LineWidth',4); 
% plot(GammaMinus_X_Cyto,GammaMinus_AVG_Cyto,'-r','LineWidth',4); 
% set(gca,'TickLength',[0.02 0.02]);
%    set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
%    box off;
%    axis([-0.25 1.25 0 150])
%    
%    ylabel('Concentration (nM)', 'fontsize', 26, 'fontweight', 'bold')
%    xlabel('Normalized time (-)', 'fontsize', 26, 'fontweight', 'bold')
%    
%       
%     hleg = legend('YPet-DnaN \gamma+','YPet-DnaN \gamma-');
% 
%    set(hleg, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%    
%      %set(gca, 'YTick', [0 10 20 30 40 50 60 70]);
%     %set(gca, 'XTick', [-0.25 0 0.25 0.5 0.75 1 1.25]);
% hold off;

%print(h3, '-dpdf', '-r600','GammaPos_GammaNeg_concentration')


%cutomeColour = ''

% h2=figure;
% 
% offsetHist = 0;
% ylabel('Probability (-)', 'fontsize', 26, 'fontweight', 'bold')
%    xlabel('Number of YPet-\beta_2(-)', 'fontsize', 26, 'fontweight', 'bold')
%    %a = annotation('textbox','Position',[0.65 0.65 0.35 0.16], 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold','FitBoxToText','on','String', sprintf('N = %d',N_CombinedSeries) );
%    axis([0 100 0 4e-2])
%    set(gca, 'XTick', [0 20 40 60 80 100]);
%    set(gca, 'YTick', [0 0.01 0.02 0.03 0.04]);
%    set(gca,'TickLength',[0.02 0.02]);
%    box off;
%    set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
%    hold on;
% bar(WT_X_hist-offsetHist/2,WT_Y_hist,'b');
% bar(gammaNegX_hist+offsetHist/2,gammaNegY_hist,'r');
% hold off;
% 
% hleg = legend('YPet-DnaN \gamma+','YPet-DnaN \gamma-');
% set(hleg, 'fontsize', 16, 'linewidth', 2, 'fontweight', 'bold');
%print(h2, '-dpdf', '-r600','GammaPos_GammaNeg_CoHIst')

%h1 = plot(axscat2D,all_perc, 'o', 'color' ,[0.4 0.4 0.4], 'MarkerSize', 2);
            
            
% plot(X_WT,WT_lo,'r-','LineWidth',2,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',5); hold on;
% plot(X_WT,WT_hi,'r-','LineWidth',2,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',5); hold on;