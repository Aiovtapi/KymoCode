clear all
close all

load('D'); %load all diffusion constants here

Nbins=100;
lb=0;
ub=2.5;
x=linspace(0,ub,Nbins);
n=4; 
A=0.5;

pdfone=@(x,D) ((n/D)^n*x.^(n-1).*exp(-n*x/D))/(factorial(n-1));
pdftwo=@(x,D1,D2) (A*(n/D1)^n*x.^(n-1).*exp(-n*x/D1))/(factorial(n-1))+(1-A)*((n/D2)^n*x.^(n-1).*exp(-n*x/D2))/(factorial(n-1));
pdftwofrac=@(x,D1,D2,B) (B*(n/D1)^n*x.^(n-1).*exp(-n*x/D1))/(factorial(n-1))+(1-B)*((n/D2)^n*x.^(n-1).*exp(-n*x/D2))/(factorial(n-1));

[phattwo,pci1]=mle(D,'pdf',pdftwo,'start',[0.01,0.01]);
[phatone,pci2]=mle(D,'pdf',pdfone,'start',0.01);
[phattwofrac,pci3]=mle(D,'pdf',pdftwofrac,'start',[0.01,0.4,0.5]);

Art=pdftwo(x,phattwo(1),phattwo(2));
Brt=pdfone(x,phatone(1));
Crt=pdftwofrac(x,phattwofrac(1),phattwofrac(2),phattwofrac(3));


Drt=pdfone(x,phattwofrac(1));
Ert=pdfone(x,phattwofrac(2));


FracLow=phattwofrac(3);
FracHigh=1-FracLow;

%% Goodness of fit
bins=linspace(lb, ub, Nbins);
Dc=histc(D,bins);


%% 
hold on
fig1=figure(1);
h=bar(bins,Dc*100/(sum(Dc)));
h.EdgeColor='k';
h.FaceColor=[.5 .5 .5];
set(gca,'fontsize',18)
axis([-0.005 2.5 0 8])
xlabel('Diffusion Coefficient (um^{2}/s)')
ylabel('Probability (-)')
% plot(x,Art/sum(Art),'r');
plot(x,phattwofrac(3)*Drt/sum(Drt)*100,'r','LineWidth',2);
plot(x,(1-phattwofrac(3))*Ert/sum(Ert)*100,'b','LineWidth',2);
title(strcat({'UvrD Diffusion PDF - Traj = '},num2str(size(D,1))));
t=annotation('textbox',[0.45 0.2 0.2 0.07],...
    'String',strcat({'D_{fast} = '},sprintf('%.3f',max(phattwofrac(1),phattwofrac(2)))),...
    'FontName','Arial');
s=annotation('textbox',[0.2 0.8 0.2 0.07],...
    'String',strcat({'D_{slow} = '},sprintf('%.3f',min(phattwofrac(1),phattwofrac(2)))),...
    'FontName','Arial');

SlowP=annotation('textbox',[0.17 0.6 0.1 0.1],...
    'String',strcat(sprintf('%2.0f',phattwofrac(3)*100),'%'),...
    'FontName','Arial',...
    'EdgeColor','w','FontSize',25,'Color','r','FontWeight','bold');
FastP=annotation('textbox',[0.3 0.2 0.1 0.1],...
    'String',strcat(sprintf('%2.0f',(1-phattwofrac(3))*100),'%'),...
    'FontName','Arial',...
    'EdgeColor','w','FontSize',25,'Color','b','FontWeight','bold');

t.FontSize = 16;
s.FontSize = 16;
t.FontWeight='bold';
s.FontWeight='bold';


