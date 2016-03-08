close all

Npoints=30;
Method=1;

DatafileName='Data';

switch Method
    
    case 1 %Combine individual frame files to a combined file
        
        Restot=zeros(Npoints,3);
        for i=1:Npoints
            n=num2str(i-1);
            Restot(i,:)=importdata(strcat('SimulationResults/KAKTEST/',n,'.mat'));
        end
        save(strcat('SimulationResults/KAKTEST/',DatafileName),'Restot');
        
    case 2 %takes two combined files and combines them to a 'Data' file.
        
        A=importdata('SimulationResults/Total.mat');
        B=importdata('SimulationResults/Total2.mat');

        Data=cat(1,A,B);
        save('SimulationResults/Data','Data');
        
    case 3 %if total data file is already there
        Restot=importdata('SimulationResults/1/Data_80discrit.mat');

end


figure(2)
hold on
scatter(Restot(:,1),Restot(:,2),100,'b','filled');
scatter(Restot(:,1),Restot(:,3),100,'r','filled');
hold off
axis([1 6 0 1.2])
legend('Sensitivity','PPV')
set(gca,'FontSize',14);
title('Sensitivity and Positive Predictive Value vs. SNR','FontSize',20)
xlabel('Signal to Noise Ratio (-)','FontSize',16,'FontWeight','bold')
ylabel('Statistical Factor (-)','FontSize',16,'FontWeight','bold')
