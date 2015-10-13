%% Load 2D Gaussian fit data and perform calculations.

clear all
close all
clc

% oriZ - Dif project
initval.basepath='/Users/rleeuw/Work/Data/20150511_OriZ_dif/Exp5/Figures/BacPics/';

%% Define variables

Ncells=4;

D=cell(Ncells+1,1);
d=cell(Ncells+1,1);
S=cell(Ncells,1);
Sd=cell(Ncells,1);
BacLife=zeros(Ncells,1); BacLifed=zeros(Ncells,1);
MainPathDnaN=strcat(initval.basepath,'DnaN2/DataMULTI/');
MainPathdif=strcat(initval.basepath,'dif2/DataMULTI/');

%% Load data
tic
for i=1:Ncells % CHANGED FOR ONE BAC 
    D{i}=load(strcat(MainPathDnaN,num2str(i),'.mat'));
    d{i}=load(strcat(MainPathdif,num2str(i),'.mat'));
end

%% Calculations and Filtering
 
 Ilb=20000;
 Ilbd=45000;

 d{4}.x(:,6)=2*d{4}.x(:,6); % Single Bac Obs
 
for i=1:Ncells
     % Filter out integrated intensities by upper bound
     % First channel
     
    Dummer=D{i}.x(:,6);
    Dummer(Dummer<Ilb)=0;
    D{i}.x(:,6)=Dummer;
    
    Dummer=D{i}.xR1(:,6);
    Dummer(Dummer<Ilb)=0;
    D{i}.xR1(:,6)=Dummer;
    
    Dummer=D{i}.xR2(:,6);
    Dummer(Dummer<Ilb)=0;
    D{i}.xR2(:,6)=Dummer;
    
    Dummer=D{i}.xR3(:,6);
    Dummer(Dummer<Ilb)=0;
    D{i}.xR3(:,6)=Dummer;
    
    % Second Channel
    
    Dummer=d{i}.x(:,6);
    Dummer(Dummer<Ilbd)=0;
    d{i}.x(:,6)=Dummer;
    
    Dummer=d{i}.xR1(:,6);
    Dummer(Dummer<Ilbd)=0;
    d{i}.xR1(:,6)=Dummer;
    
    Dummer=d{i}.xR2(:,6);
    Dummer(Dummer<Ilbd)=0;
    d{i}.xR2(:,6)=Dummer;
    
    Dummer=d{i}.xR3(:,6);
    Dummer(Dummer<Ilbd)=0;
    d{i}.xR3(:,6)=Dummer;
    
end

n=Ncells+1;

%% Spot tracking calculations TUS -- switch positions between spots given criterium -- switch intensities
% exchange the spot positions with ones closest to previous values

for i=1:Ncells
    for j=2:length(D{i}.x(:,1));
        
        pixelsbetwspots=4;
        CritValue=pixelsbetwspots/D{i}.XSize(j);

        Diffvectorx=[abs(D{i}.XNorm(j,2)-D{i}.XNorm(j-1,2)) ...
            abs(D{i}.XNormR1(j,2)-D{i}.XNorm(j-1,2)) ...
            abs(D{i}.XNormR2(j,2)-D{i}.XNorm(j-1,2)) ...
            abs(D{i}.XNormR3(j,2)-D{i}.XNorm(j-1,2))];

        if D{i}.XNormR1(j,2)==0
            Diffvectorx(2)=2;
        end
        if D{i}.XNormR2(j,2)==0
            Diffvectorx(3)=2;
        end
        if D{i}.XNormR3(j,2)==0
            Diffvectorx(4)=2;
        end    
        
idx=find(Diffvectorx<=CritValue);
% 
% figure(1)
% hold on 
% plot(D{i}.XNorm(j,2),0.5,'ob','MarkerSize',16,'LineWidth',2);
% plot(D{i}.XNormR1(j,2),0.5,'or','MarkerSize',14,'LineWidth',2);
% plot(D{i}.XNormR2(j,2),0.5,'ok','MarkerSize',12,'LineWidth',2);
% plot(D{i}.XNormR3(j,2),0.5,'og','MarkerSize',10,'LineWidth',2);
% plot(D{i}.XNorm(j-1,2),0,'ob','MarkerSize',16,'LineWidth',2);
% plot(D{i}.XNormR1(j-1,2),0,'or','MarkerSize',14,'LineWidth',2);
% plot(D{i}.XNormR2(j-1,2),0,'ok','MarkerSize',12,'LineWidth',2);
% plot(D{i}.XNormR3(j-1,2),0,'og','MarkerSize',10,'LineWidth',2);
% axis([0 1 0.9 1.2])

% If spots merge, total intensity goes to 'bigger' spot.

    if sum(idx)==2;
    DUM=D{i}.XNorm(j,2);
    
    D{i}.XNorm(j,2)=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=D{i}.XNorm(j,2);
    
    D{i}.XNorm(j,2)=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(D{i}.XNormR1(j,2)+D{i}.XNorm(j,2))/2;
    
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=D{i}.XNorm(j,2);
    D{i}.XNorm(j,2)=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(D{i}.XNorm(j,2)+D{i}.XNormR2(j,2))/2;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(D{i}.XNorm(j,2)+D{i}.XNormR3(j,2))/2;
    
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;

    elseif sum(idx)==5 && idx(1)==2;
    DUM=D{i}.XNorm(j,2);
    DUM2=(D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM2;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==2;
    DUM=D{i}.XNorm(j,2);
    DUM2=(D{i}.XNormR1(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM2;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==1;
    DUM=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/3;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==7;
    DUM=D{i}.XNorm(j,2);
    DUM2=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNorm(j,2)=DUM2;
    D{i}.XNormR3(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==8;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNorm(j,2))/3;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=D{i}.XNorm(j,2);
    DUM2=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2))/3;
    D{i}.XNorm(j,2)=DUM2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==10;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2)+D{i}.XNorm(j,2))/4;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    end

% plot(D{i}.XNorm(j,2),1,'ob','MarkerSize',16,'LineWidth',2);
% plot(D{i}.XNormR1(j,2),1,'or','MarkerSize',14,'LineWidth',2);
% plot(D{i}.XNormR2(j,2),1,'ok','MarkerSize',12,'LineWidth',2);
% plot(D{i}.XNormR3(j,2),1,'og','MarkerSize',10,'LineWidth',2);

% Second Spot ---------------------------------------------------------------------

% for this one I send R2 to the other (distant) spot, or R3, if normal x not there.

DiffvectorxR1=[abs(D{i}.XNormR1(j,2)-D{i}.XNormR1(j-1,2)) ...
            abs(D{i}.XNorm(j,2)-D{i}.XNormR1(j-1,2)) ...
            abs(D{i}.XNormR2(j,2)-D{i}.XNormR1(j-1,2)) ...
            abs(D{i}.XNormR3(j,2)-D{i}.XNormR1(j-1,2))];

        
        if D{i}.XNorm(j,2)==0
            DiffvectorxR1(2)=2;
        end
        if D{i}.XNormR2(j,2)==0
            DiffvectorxR1(3)=2;
        end
        if D{i}.XNormR3(j,2)==0
            DiffvectorxR1(4)=2;
        end   
        
idx=find(DiffvectorxR1<=CritValue);

    if sum(idx)==2;
    DUM=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=D{i}.XNorm(j,2);
    D{i}.XNorm(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(D{i}.XNormR1(j,2)+D{i}.XNorm(j,2))/2;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(D{i}.XNormR1(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==2;
    DUM=D{i}.XNormR1(j,2);
    DUM2=(D{i}.XNorm(j,2)+D{i}.XNormR2(j,2))/2;
    D{i}.XNormR1(j,2)=DUM2;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==2;
    DUM=D{i}.XNormR1(j,2);
    DUM2=(D{i}.XNorm(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR1(j,2)=DUM2;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==1;
    DUM=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/3;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==7;
    DUM=D{i}.XNormR1(j,2);
    DUM2=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR1(j,2)=DUM2;
    D{i}.XNormR3(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==8;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2))/3;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=D{i}.XNormR1(j,2);
    DUM2=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNorm(j,2))/3;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM2;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==10;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2)+D{i}.XNorm(j,2))/4;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    end
    
% plot(D{i}.XNorm(j,2),2,'ob','MarkerSize',16,'LineWidth',2);
% plot(D{i}.XNormR1(j,2),2,'or','MarkerSize',14,'LineWidth',2);
% plot(D{i}.XNormR2(j,2),2,'ok','MarkerSize',12,'LineWidth',2);
% plot(D{i}.XNormR3(j,2),2,'og','MarkerSize',10,'LineWidth',2);

% THIRD SPOT -------------------------------------------------


DiffvectorxR2=[abs(D{i}.XNormR2(j,2)-D{i}.XNormR2(j-1,2)) ... 
            abs(D{i}.XNorm(j,2)-D{i}.XNormR2(j-1,2)) ...
            abs(D{i}.XNormR1(j,2)-D{i}.XNormR2(j-1,2)) ...
            abs(D{i}.XNormR3(j,2)-D{i}.XNormR2(j-1,2))];
        
        if D{i}.XNormR2(j,2)==0
            DiffvectorxR2(1:4)=2;
            D{i}.XNormR2(j,2)=D{i}.XNormR2(j-1,2);
        end
        if D{i}.XNorm(j,2)==0
            DiffvectorxR2(2)=2;
        end
        if D{i}.XNormR1(j,2)==0
            DiffvectorxR2(3)=2;
        end
        if D{i}.XNormR3(j,2)==0
            DiffvectorxR2(4)=2;
        end   

idx=find(DiffvectorxR2<=CritValue);

    if sum(idx)==2;
    DUM=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=D{i}.XNorm(j,2);
    D{i}.XNorm(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNorm(j,2))/2;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==2;
    DUM=D{i}.XNormR2(j,2);
    DUM2=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM2;
    
    elseif sum(idx)==6 && idx(1)==2;
    DUM=D{i}.XNormR2(j,2);
    DUM2=(D{i}.XNorm(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR2(j,2)=DUM2;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==1;
    DUM=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/3;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==7;
    DUM=D{i}.XNormR2(j,2);
    DUM2=(D{i}.XNormR1(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM2;
    
    elseif sum(idx)==8;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2))/3;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=D{i}.XNormR2(j,2);
    DUM2=(D{i}.XNormR1(j,2)+D{i}.XNormR3(j,2)+D{i}.XNorm(j,2))/3;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM2;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==10;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2)+D{i}.XNorm(j,2))/4;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    end
    
% plot(D{i}.XNorm(j,2),3,'ob','MarkerSize',16,'LineWidth',2);
% plot(D{i}.XNormR1(j,2),3,'or','MarkerSize',14,'LineWidth',2);
% plot(D{i}.XNormR2(j,2),3,'ok','MarkerSize',12,'LineWidth',2);
% plot(D{i}.XNormR3(j,2),3,'og','MarkerSize',10,'LineWidth',2);

% FOURTH SPOT --------------------------------------------------------

DiffvectorxR3=[abs(D{i}.XNormR3(j,2)-D{i}.XNormR3(j-1,2))...
            abs(D{i}.XNorm(j,2)-D{i}.XNormR3(j-1,2)) ...
            abs(D{i}.XNormR1(j,2)-D{i}.XNormR3(j-1,2)) ...
            abs(D{i}.XNormR2(j,2)-D{i}.XNormR3(j-1,2))];
        
        if D{i}.XNormR3(j,2)==0
            DiffvectorxR3(1:4)=2;
            D{i}.XNormR3(j,2)=D{i}.XNormR3(j-1,2);
        end
        if D{i}.XNorm(j,2)==0
            DiffvectorxR3(2)=2;
        end
        if D{i}.XNormR1(j,2)==0
            DiffvectorxR3(3)=2;
        end
        if D{i}.XNormR2(j,2)==0
            DiffvectorxR3(4)=2;
        end   

idx=find(DiffvectorxR3<=CritValue);

    if sum(idx)==2;
    DUM=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=D{i}.XNorm(j,2);
    D{i}.XNorm(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=D{i}.XNormR1(j,2);
    D{i}.XNormR1(j,2)=DUM; 
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(D{i}.XNormR3(j,2)+D{i}.XNorm(j,2))/2;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=D{i}.XNormR3(j,2);
    D{i}.XNormR3(j,2)=D{i}.XNormR2(j,2);
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(D{i}.XNormR1(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2))/2;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==2;
    DUM=D{i}.XNormR3(j,2);
    DUM2=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM2;
    
    elseif sum(idx)==6 && idx(1)==1;
    DUM=(D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/3;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==2;
    DUM=D{i}.XNormR3(j,2);
    DUM2=(D{i}.XNorm(j,2)+D{i}.XNormR2(j,2))/2;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM2;
    
    elseif sum(idx)==7;
    DUM=D{i}.XNormR3(j,2);
    DUM2=(D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/2;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM2;
    D{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==8;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2))/3;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=D{i}.XNormR3(j,2);
    DUM2=(D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2)+D{i}.XNorm(j,2))/3;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM2;
    
    elseif sum(idx)==10;
    DUM=(D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2)+D{i}.XNormR1(j,2)+D{i}.XNorm(j,2))/4;
    D{i}.XNorm(j,2)=DUM;
    D{i}.XNormR1(j,2)=DUM;
    D{i}.XNormR2(j,2)=DUM;
    D{i}.XNormR3(j,2)=DUM;
    end
    
% plot(D{i}.XNorm(j,2),4,'ob','MarkerSize',16,'LineWidth',2);
% plot(D{i}.XNormR1(j,2),4,'or','MarkerSize',14,'LineWidth',2);
% plot(D{i}.XNormR2(j,2),4,'ok','MarkerSize',12,'LineWidth',2);
% plot(D{i}.XNormR3(j,2),4,'og','MarkerSize',10,'LineWidth',2);
% hold off
% axis([0 1 0 4.5])
% legend('Brightest','R1','R2','R3');
% set(gca,'FontSize',16);
% xlabel('Normalized Cell Position (-)');
% ylabel('Algorithm Chronology (-)');
    end
end


% SECOND CHANNEL ------------------------------------------
for i=1:Ncells
    for j=2:length(d{i}.x(:,1));
        
        pixelsbetwspots=6;
        CritValue=pixelsbetwspots/d{i}.XSize(j);
        
        Diffvectorx=[abs(d{i}.XNorm(j,2)-d{i}.XNorm(j-1,2)) ...
            abs(d{i}.XNormR1(j,2)-d{i}.XNorm(j-1,2)) ...
            abs(d{i}.XNormR2(j,2)-d{i}.XNorm(j-1,2)) ...
            abs(d{i}.XNormR3(j,2)-d{i}.XNorm(j-1,2))];
       
% Make sure that non existing points are not considered.

        if d{i}.XNormR1(j,2)==0
            Diffvectorx(2)=2;
        end
        if d{i}.XNormR2(j,2)==0
            Diffvectorx(3)=2;
        end
        if d{i}.XNormR3(j,2)==0
            Diffvectorx(4)=2;
        end            
           
idx=find(Diffvectorx<=CritValue);


% If summed index is 2, current (j) XNorm is switched with XNormR1 and vv.
% If summed index is 3, current (j) XNorm is switched with XNormR2 , 
% or if index is 1+2 XNorm and XNormR1 are changed to the COM of XNorm and
% XNormR1. etc.
% If summed index is such that XNorm is COM of two other (e.g. R1 & R2)
% spots, then one of other spots goes to XNorm position. I always choose
% R1, or if not R1 is not considered, R2 to go to the other distant spot.
j
    if sum(idx)==2;
    DUM=d{i}.XNorm(j,2);
    d{i}.XNorm(j,2)=d{i}.XNormR1(j,2);
    d{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=d{i}.XNorm(j,2);
    d{i}.XNorm(j,2)=d{i}.XNormR2(j,2);
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(d{i}.XNormR1(j,2)+d{i}.XNorm(j,2))/2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=d{i}.XNorm(j,2);
    d{i}.XNorm(j,2)=d{i}.XNormR3(j,2);
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(d{i}.XNorm(j,2)+d{i}.XNormR2(j,2))/2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(d{i}.XNorm(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==2;
    DUM=d{i}.XNorm(j,2);
    DUM2=(d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM2;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==1;
    DUM=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/3;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;

    elseif sum(idx)==6 && idx(1)==2;
    DUM=d{i}.XNorm(j,2);
    DUM2=(d{i}.XNormR1(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM2;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==7;
    DUM=d{i}.XNorm(j,2);
    DUM2=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNorm(j,2)=DUM2;
    d{i}.XNormR3(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==8;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNorm(j,2))/3;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=d{i}.XNorm(j,2);
    DUM2=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2))/3;
    d{i}.XNorm(j,2)=DUM2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==10;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2)+d{i}.XNorm(j,2))/4;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    end

% Second Spot ---------------------------------------------------------------------

% for this one I send R2 to the other (distant) spot, or R3, if normal x not there.

DiffvectorxR1=[abs(d{i}.XNormR1(j,2)-d{i}.XNormR1(j-1,2)) ...
            abs(d{i}.XNorm(j,2)-d{i}.XNormR1(j-1,2)) ...
            abs(d{i}.XNormR2(j,2)-d{i}.XNormR1(j-1,2)) ...
            abs(d{i}.XNormR3(j,2)-d{i}.XNormR1(j-1,2))];

% Make sure that non existing points are not considered.

        if d{i}.XNorm(j,2)==0
            DiffvectorxR1(2)=2;
        end
        if d{i}.XNormR2(j,2)==0
            DiffvectorxR1(3)=2;
        end
        if d{i}.XNormR3(j,2)==0
            DiffvectorxR1(4)=2;
        end    
        
idx=find(DiffvectorxR1<=CritValue);

    if sum(idx)==2;
    DUM=d{i}.XNormR1(j,2);
    d{i}.XNormR1(j,2)=d{i}.XNorm(j,2);
    d{i}.XNorm(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=d{i}.XNormR1(j,2);
    d{i}.XNormR1(j,2)=d{i}.XNormR2(j,2);
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(d{i}.XNormR1(j,2)+d{i}.XNorm(j,2))/2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=d{i}.XNormR1(j,2);
    d{i}.XNormR1(j,2)=d{i}.XNormR3(j,2);
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(d{i}.XNormR1(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==2;
    DUM=d{i}.XNormR1(j,2);
    DUM2=(d{i}.XNorm(j,2)+d{i}.XNormR2(j,2))/2;
    d{i}.XNormR1(j,2)=DUM2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;

    elseif sum(idx)==6 && idx(1)==1;
    DUM=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/3;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==2;
    DUM=d{i}.XNormR1(j,2);
    DUM2=(d{i}.XNorm(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR1(j,2)=DUM2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==7;
    DUM=d{i}.XNormR1(j,2);
    DUM2=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR1(j,2)=DUM2;
    d{i}.XNormR3(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==8;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2))/3;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=d{i}.XNormR1(j,2);
    DUM2=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNorm(j,2))/3;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM2;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==10;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2)+d{i}.XNorm(j,2))/4;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    end

    
    
% THIRD SPOT -------------------------------------------------


DiffvectorxR2=[abs(d{i}.XNormR2(j,2)-d{i}.XNormR2(j-1,2)) ... 
            abs(d{i}.XNorm(j,2)-d{i}.XNormR2(j-1,2)) ...
            abs(d{i}.XNormR1(j,2)-d{i}.XNormR2(j-1,2)) ...
            abs(d{i}.XNormR3(j,2)-d{i}.XNormR2(j-1,2))];
        
% Make sure that non existing points are not considered.

        if d{i}.XNormR2(j,2)==0
            DiffvectorxR2(1:4)=2;
            d{i}.XNormR2(j,2)=d{i}.XNormR2(j-1,2);
        end
        if d{i}.XNorm(j,2)==0
            DiffvectorxR2(2)=2;
        end
        if d{i}.XNormR1(j,2)==0
            DiffvectorxR2(3)=2;
        end
        if d{i}.XNormR3(j,2)==0
            DiffvectorxR2(4)=2;
        end  

idx=find(DiffvectorxR2<=CritValue);

if  sum(idx)==2;
    DUM=d{i}.XNormR2(j,2);
    d{i}.XNormR2(j,2)=d{i}.XNorm(j,2);
    d{i}.XNorm(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=d{i}.XNormR2(j,2);
    d{i}.XNormR2(j,2)=d{i}.XNormR1(j,2);
    d{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNorm(j,2))/2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=d{i}.XNormR2(j,2);
    d{i}.XNormR2(j,2)=d{i}.XNormR3(j,2);
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==2;
    DUM=d{i}.XNormR2(j,2);
    DUM2=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM2;
    d{i}.XNormR2(j,2)=DUM2;
    
    elseif sum(idx)==6 && idx(1)==1;
    DUM=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/3;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==2;
    DUM=d{i}.XNormR2(j,2);
    DUM2=(d{i}.XNorm(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR2(j,2)=DUM2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==7;
    DUM=d{i}.XNormR2(j,2);
    DUM2=(d{i}.XNormR1(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM2;
    
    elseif sum(idx)==8;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2))/3;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=d{i}.XNormR2(j,2);
    DUM2=(d{i}.XNormR1(j,2)+d{i}.XNormR3(j,2)+d{i}.XNorm(j,2))/3;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM2;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==10;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2)+d{i}.XNorm(j,2))/4;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
end

% FOURTH SPOT --------------------------------------------------------

DiffvectorxR3=[abs(d{i}.XNormR3(j,2)-d{i}.XNormR3(j-1,2))...
            abs(d{i}.XNorm(j,2)-d{i}.XNormR3(j-1,2)) ...
            abs(d{i}.XNormR1(j,2)-d{i}.XNormR3(j-1,2)) ...
            abs(d{i}.XNormR2(j,2)-d{i}.XNormR3(j-1,2))];

% Make sure that non existing points are not considered.

        if d{i}.XNormR3(j,2)==0
            DiffvectorxR3(1:4)=2;
            d{i}.XNormR3(j,2)=d{i}.XNormR3(j-1,2);
        end
        if d{i}.XNorm(j,2)==0
            DiffvectorxR3(2)=2;
        end
        if d{i}.XNormR1(j,2)==0
            DiffvectorxR3(3)=2;
        end
        if d{i}.XNormR2(j,2)==0
            DiffvectorxR3(4)=2;
        end  
        
        
idx=find(DiffvectorxR3<=CritValue);

    if sum(idx)==2;
    DUM=d{i}.XNormR3(j,2);
    d{i}.XNormR3(j,2)=d{i}.XNorm(j,2);
    d{i}.XNorm(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==3;
    DUM=d{i}.XNormR3(j,2);
    d{i}.XNormR3(j,2)=d{i}.XNormR1(j,2);
    d{i}.XNormR1(j,2)=DUM;
    
    elseif sum(idx)==3 && idx(1)==1;
    DUM=(d{i}.XNormR3(j,2)+d{i}.XNorm(j,2))/2;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==4;
    DUM=d{i}.XNormR3(j,2);
    d{i}.XNormR3(j,2)=d{i}.XNormR2(j,2);
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==4 && idx(1)==1;
    DUM=(d{i}.XNormR1(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==1;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2))/2;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==5 && idx(1)==2;
    DUM=d{i}.XNormR3(j,2);
    DUM2=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM2;
    
    elseif sum(idx)==6 && idx(1)==1;
    DUM=(d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/3;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==6 && idx(1)==2;
    DUM=d{i}.XNormR3(j,2);
    DUM2=(d{i}.XNorm(j,2)+d{i}.XNormR2(j,2))/2;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM2;
    
    elseif sum(idx)==7;
    DUM=d{i}.XNormR3(j,2);
    DUM2=(d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/2;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM2;
    d{i}.XNormR2(j,2)=DUM;
    
    elseif sum(idx)==8;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2))/3;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    
    elseif sum(idx)==9;
    DUM=d{i}.XNormR3(j,2);
    DUM2=(d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2)+d{i}.XNorm(j,2))/3;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM2;
    
    elseif sum(idx)==10;
    DUM=(d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2)+d{i}.XNormR1(j,2)+d{i}.XNorm(j,2))/4;
    d{i}.XNorm(j,2)=DUM;
    d{i}.XNormR1(j,2)=DUM;
    d{i}.XNormR2(j,2)=DUM;
    d{i}.XNormR3(j,2)=DUM;
    end
    
    end
end

%% Also when spots are close together, position and II merge.

% for i=1:Ncells
%     for j=1:length(D{i}.XNorm(:,2));      
%         
%         pixelsbetwspots=7;
%         CritValueMerge=pixelsbetwspots/D{i}.XSize(j);
% 
%         Diffvectorx=[abs(D{i}.XNormR1(j,2)-D{i}.XNorm(j,2)) ...
%             abs(D{i}.XNormR2(j,2)-D{i}.XNorm(j,2)) ...
%             abs(D{i}.XNormR3(j,2)-D{i}.XNorm(j,2))];
%         
% idx=find(Diffvectorx<=CritValueMerge);
% 
% if sum(idx)==1
%     DUM=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2))/2;
%     D{i}.XNormR1(j,2)=DUM;
%     D{i}.XNorm(j,2)=DUM;
%     D{i}.x(j,6)=D{i}.x(j,6)+D{i}.xR1(j,6);
%     D{i}.xR1(j,6)=0;
% elseif sum(idx)==2
%     DUM=(D{i}.XNorm(j,2)+D{i}.XNormR2(j,2))/2;
%     D{i}.XNormR2(j,2)=DUM;
%     D{i}.XNorm(j,2)=DUM;
%     D{i}.x(j,6)=D{i}.x(j,6)+D{i}.xR2(j,6);
%     D{i}.xR2(j,6)=0;    
% elseif sum(idx)==3 && idx(1)==3
%     DUM=(D{i}.XNorm(j,2)+D{i}.XNormR3(j,2))/2;
%     D{i}.XNormR3(j,2)=DUM;
%     D{i}.XNorm(j,2)=DUM;
%     D{i}.x(j,6)=D{i}.x(j,6)+D{i}.xR3(j,6);
%     D{i}.xR3(j,6)=0;  
% elseif sum(idx)==3 && idx(1)==1
%     DUM=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2))/3;
%     D{i}.XNormR1(j,2)=DUM;
%     D{i}.XNormR2(j,2)=DUM;
%     D{i}.XNorm(j,2)=DUM;
%     D{i}.x(j,6)=D{i}.x(j,6)+D{i}.xR2(j,6)+D{i}.xR1(j,6);
%     D{i}.xR2(j,6)=0; D{i}.xR1(j,6)=0;
% elseif sum(idx)==4 
%     DUM=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2)+D{i}.XNormR3(j,2))/3;
%     D{i}.XNormR1(j,2)=DUM;
%     D{i}.XNormR3(j,2)=DUM;
%     D{i}.XNorm(j,2)=DUM;
%     D{i}.x(j,6)=D{i}.x(j,6)+D{i}.xR1(j,6)+D{i}.xR3(j,6);
%     D{i}.xR3(j,6)=0; D{i}.xR1(j,6)=0;
% elseif sum(idx)==5
%     DUM=(D{i}.XNorm(j,2)+D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2))/3;
%     D{i}.XNormR2(j,2)=DUM;
%     D{i}.XNormR3(j,2)=DUM;
%     D{i}.XNorm(j,2)=DUM;
%     D{i}.x(j,6)=D{i}.x(j,6)+D{i}.xR2(j,6)+D{i}.xR3(j,6);
%     D{i}.xR3(j,6)=0; D{i}.xR2(j,6)=0;
%  elseif sum(idx)==6
%     DUM=(D{i}.XNorm(j,2)+D{i}.XNormR1(j,2)+D{i}.XNormR2(j,2)+D{i}.XNormR3(j,2))/4;
%     D{i}.XNormR1(j,2)=DUM;
%     D{i}.XNormR2(j,2)=DUM;
%     D{i}.XNormR3(j,2)=DUM;
%     D{i}.XNorm(j,2)=DUM;
%     D{i}.x(j,6)=D{i}.x(j,6)+D{i}.xR1(j,6)+D{i}.xR2(j,6)+D{i}.xR3(j,6);
%     D{i}.xR3(j,6)=0; D{i}.xR2(j,6)=0; D{i}.xR1(j,6)=0;
% end
%     end   
%     
%     for j=1:length(d{i}.XNorm(:,2));      
% 
%                 
%         pixelsbetwspots=4;
%         CritValueMerge=pixelsbetwspots/D{i}.XSize(j);
% 
%         Diffvectorx=[abs(d{i}.XNormR1(j,2)-d{i}.XNorm(j,2)) ...
%             abs(d{i}.XNormR2(j,2)-d{i}.XNorm(j,2)) ...
%             abs(d{i}.XNormR3(j,2)-d{i}.XNorm(j,2))];
%         
% idx=find(Diffvectorx<=CritValueMerge);
% 
% if sum(idx)==1
%     DUM=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2))/2;
%     d{i}.XNormR1(j,2)=DUM;
%     d{i}.XNorm(j,2)=DUM;
%     d{i}.x(j,6)=d{i}.x(j,6)+d{i}.xR1(j,6);
%     d{i}.xR1(j,6)=0;
% elseif sum(idx)==2
%     DUM=(d{i}.XNorm(j,2)+d{i}.XNormR2(j,2))/2;
%     d{i}.XNormR2(j,2)=DUM;
%     d{i}.XNorm(j,2)=DUM;
%     d{i}.x(j,6)=d{i}.x(j,6)+d{i}.xR2(j,6);
%     d{i}.xR2(j,6)=0;    
% elseif sum(idx)==3 && idx(1)==3
%     DUM=(d{i}.XNorm(j,2)+d{i}.XNormR3(j,2))/2;
%     d{i}.XNormR3(j,2)=DUM;
%     d{i}.XNorm(j,2)=DUM;
%     d{i}.x(j,6)=d{i}.x(j,6)+d{i}.xR3(j,6);
%     d{i}.xR3(j,6)=0;  
% elseif sum(idx)==3 && idx(1)==1
%     DUM=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2))/3;
%     d{i}.XNormR1(j,2)=DUM;
%     d{i}.XNormR2(j,2)=DUM;
%     d{i}.XNorm(j,2)=DUM;
%     d{i}.x(j,6)=d{i}.x(j,6)+d{i}.xR2(j,6)+d{i}.xR1(j,6);
%     d{i}.xR2(j,6)=0; d{i}.xR1(j,6)=0;
% elseif sum(idx)==4 
%     DUM=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2)+d{i}.XNormR3(j,2))/3;
%     d{i}.XNormR1(j,2)=DUM;
%     d{i}.XNormR3(j,2)=DUM;
%     d{i}.XNorm(j,2)=DUM;
%     d{i}.x(j,6)=d{i}.x(j,6)+d{i}.xR1(j,6)+d{i}.xR3(j,6);
%     d{i}.xR3(j,6)=0; d{i}.xR1(j,6)=0;
% elseif sum(idx)==5
%     DUM=(d{i}.XNorm(j,2)+d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2))/3;
%     d{i}.XNormR2(j,2)=DUM;
%     d{i}.XNormR3(j,2)=DUM;
%     d{i}.XNorm(j,2)=DUM;
%     d{i}.x(j,6)=d{i}.x(j,6)+d{i}.xR2(j,6)+d{i}.xR3(j,6);
%     d{i}.xR3(j,6)=0; d{i}.xR2(j,6)=0;
% elseif sum(idx)==6
%     DUM=(d{i}.XNorm(j,2)+d{i}.XNormR1(j,2)+d{i}.XNormR2(j,2)+d{i}.XNormR3(j,2))/4;
%     d{i}.XNormR1(j,2)=DUM;
%     d{i}.XNormR2(j,2)=DUM;
%     d{i}.XNormR3(j,2)=DUM;
%     d{i}.XNorm(j,2)=DUM;
%     d{i}.x(j,6)=d{i}.x(j,6)+d{i}.xR1(j,6)+d{i}.xR2(j,6)+d{i}.xR3(j,6);
%     d{i}.xR3(j,6)=0; d{i}.xR2(j,6)=0; d{i}.xR1(j,6)=0;
% end
% 
%     end
% end

%% Calcs

% Initiate Spot Intensity Cells

D{n}.I=D{1}.x(:,1);
D{n}.IR1=D{1}.xR1(:,1);
D{n}.IR2=D{1}.xR2(:,1);
D{n}.IR3=D{1}.xR3(:,1);

d{n}.I=d{1}.x(:,1);
d{n}.IR1=d{1}.xR1(:,1);
d{n}.IR2=d{1}.xR2(:,1);
d{n}.IR3=d{1}.xR3(:,1);

% Initiate Cell Intensity Cells

D{n}.FCI=D{1}.x(:,7);

d{n}.FCI=d{1}.x(:,7);

% Initiate X Position Cells

D{n}.X=D{1}.XNorm(:,2);
D{n}.XR1=D{1}.XNormR1(:,2);
D{n}.XR2=D{1}.XNormR2(:,2);
D{n}.XR3=D{1}.XNormR3(:,2);

d{n}.X=d{1}.XNorm(:,2);
d{n}.XR1=d{1}.XNormR1(:,2);
d{n}.XR2=d{1}.XNormR2(:,2);
d{n}.XR3=d{1}.XNormR3(:,2);

% Initiate Y Position Cells

D{n}.Y=D{1}.XNorm(:,4); 
D{n}.YR1=D{1}.XNormR1(:,4);
D{n}.YR2=D{1}.XNormR2(:,4);
D{n}.YR3=D{1}.XNormR3(:,4);

d{n}.Y=d{1}.XNorm(:,4);
d{n}.YR1=d{1}.XNormR1(:,4);
d{n}.YR2=d{1}.XNormR2(:,4);
d{n}.YR3=d{1}.XNormR3(:,4);

% Initiate Integrated Intensity Cells

D{n}.IntI=D{1}.x(:,6);
D{n}.IntIR1=D{1}.xR1(:,6);
D{n}.IntIR2=D{1}.xR2(:,6);
D{n}.IntIR3=D{1}.xR3(:,6);

d{n}.IntI=d{1}.x(:,6);
d{n}.IntIR1=d{1}.xR1(:,6);
d{n}.IntIR2=d{1}.xR2(:,6);
d{n}.IntIR3=d{1}.xR3(:,6);

for i=2:Ncells
D{n}.I=cat(1,D{n}.I,D{i}.x(:,1));
D{n}.IR1=cat(1,D{n}.IR1,D{i}.xR1(:,1));
D{n}.IR2=cat(1,D{n}.IR2,D{i}.xR2(:,1));
D{n}.IR3=cat(1,D{n}.IR3,D{i}.xR3(:,1));

D{n}.FCI=cat(1,D{n}.FCI,D{i}.x(:,7));

d{n}.FCI=cat(1,d{n}.FCI,d{i}.x(:,7));

d{n}.I=cat(1,d{n}.I,d{i}.x(:,1));
d{n}.IR1=cat(1,d{n}.IR1,d{i}.xR1(:,1));
d{n}.IR2=cat(1,d{n}.IR2,d{i}.xR2(:,1));
d{n}.IR3=cat(1,d{n}.IR3,d{i}.xR3(:,1));

D{n}.X=cat(1,D{n}.X,D{i}.XNorm(:,2));
D{n}.XR1=cat(1,D{n}.XR1,D{i}.XNormR1(:,2));
D{n}.XR2=cat(1,D{n}.XR2,D{i}.XNormR2(:,2));
D{n}.XR3=cat(1,D{n}.XR3,D{i}.XNormR3(:,2));

d{n}.X=cat(1,d{n}.X,d{i}.XNorm(:,2));
d{n}.XR1=cat(1,d{n}.XR1,d{i}.XNormR1(:,2));
d{n}.XR2=cat(1,d{n}.XR2,d{i}.XNormR2(:,2));
d{n}.XR3=cat(1,d{n}.XR3,d{i}.XNormR3(:,2));

D{n}.Y=cat(1,D{n}.Y,D{i}.XNorm(:,4));
D{n}.YR1=cat(1,D{n}.YR1,D{i}.XNormR1(:,4));
D{n}.YR2=cat(1,D{n}.YR2,D{i}.XNormR2(:,4));
D{n}.YR3=cat(1,D{n}.YR3,D{i}.XNormR3(:,4));

d{n}.Y=cat(1,d{n}.Y,d{i}.XNorm(:,4));
d{n}.YR1=cat(1,d{n}.YR1,d{i}.XNormR1(:,4));
d{n}.YR2=cat(1,d{n}.YR2,d{i}.XNormR2(:,4));
d{n}.YR3=cat(1,d{n}.YR3,d{i}.XNormR3(:,4));

D{n}.IntI=cat(1,D{n}.IntI,D{i}.x(:,6));
D{n}.IntIR1=cat(1,D{n}.IntIR1,D{i}.xR1(:,6));
D{n}.IntIR2=cat(1,D{n}.IntIR2,D{i}.xR2(:,6));
D{n}.IntIR3=cat(1,D{n}.IntIR3,D{i}.xR3(:,6));

d{n}.IntI=cat(1,d{n}.IntI,d{i}.x(:,6));
d{n}.IntIR1=cat(1,d{n}.IntIR1,d{i}.xR1(:,6));
d{n}.IntIR2=cat(1,d{n}.IntIR2,d{i}.xR2(:,6));
d{n}.IntIR3=cat(1,d{n}.IntIR3,d{i}.xR3(:,6));
end

%% Some Calculations
TotCellsStr=sprintf('Ncells = %d',Ncells);

MeanIntD=mean(mean(D{n}.IntI));
MeanIntDR1=mean(mean(D{n}.IntIR1));
MeanIntDR2=mean(mean(D{n}.IntIR2));
MeanIntDR3=mean(mean(D{n}.IntIR3));

MeanIntFCD=mean(D{n}.FCI);

StdIntD=std(std(D{n}.IntI));
StdIntDR1=std(std(D{n}.IntIR1));
StdIntDR2=std(std(D{n}.IntIR2));
StdIntDR3=std(std(D{n}.IntIR3));

StdFCIntD=std(D{n}.FCI);

MeanIntStrD=sprintf('Mean = %g',MeanIntD);
MeanIntStrDR1=sprintf('Mean = %g',MeanIntDR1);
MeanIntStrDR2=sprintf('Mean = %g',MeanIntDR2);
MeanIntStrDR3=sprintf('Mean = %g',MeanIntDR3);

StdIntStrD=sprintf('std = %g',StdIntD);
StdIntStrDR1=sprintf('std = %g',StdIntDR1);
StdIntStrDR2=sprintf('std = %g',StdIntDR2);
StdIntStrDR3=sprintf('std = %g',StdIntDR3);

MeanIntd=mean(d{n}.IntI);
MeanIntdR1=mean(d{n}.IntIR1);
MeanIntdR2=mean(d{n}.IntIR2);
MeanIntdR3=mean(d{n}.IntIR3);

MeanIntFCd=mean(d{n}.FCI);

StdIntd=std(d{n}.IntI);
StdIntdR1=std(d{n}.IntIR1);
StdIntdR2=std(d{n}.IntIR2);
StdIntdR3=std(d{n}.IntIR3);

StdFCIntd=std(d{n}.FCI);

MeanIntStrd=sprintf('Mean = %g',MeanIntd);
MeanIntStrdR1=sprintf('Mean = %g',MeanIntdR1);
MeanIntStrdR2=sprintf('Mean = %g',MeanIntdR2);
MeanIntStrdR3=sprintf('Mean = %g',MeanIntdR3);

StdIntStrd=sprintf('std = %g',StdIntd);
StdIntStrdR1=sprintf('std = %g',StdIntdR1);
StdIntStrdR2=sprintf('std = %g',StdIntdR2);
StdIntStrdR3=sprintf('std = %g',StdIntdR3);

for i=1:Ncells
    BacLife(i)=length(D{i}.x(:,1));
    BacLifed(i)=length(d{i}.x(:,1));
end

MeanBacLife=mean(BacLife);
MeanBacLifed=mean(BacLifed);

MaxBacLife=max(BacLife);
MaxBacLifed=max(BacLifed);

% Define matrices storing respective parameter (e.g. intensity)
% for different frames (rows) for all cells (columns)

Ki=zeros(MeanBacLife,Ncells);
KiFC=zeros(MeanBacLife,Ncells);
KiR1=zeros(MaxBacLife,Ncells);
KiR2=zeros(MaxBacLife,Ncells);
KiR3=zeros(MaxBacLife,Ncells);

Kx=zeros(MeanBacLife,Ncells);
KxR1=zeros(MeanBacLife,Ncells);
KxR2=zeros(MeanBacLife,Ncells);
KxR3=zeros(MeanBacLife,Ncells);

Ky=zeros(MeanBacLife,Ncells);
KyR1=zeros(MeanBacLife,Ncells);
KyR2=zeros(MeanBacLife,Ncells);
KyR3=zeros(MeanBacLife,Ncells);

Kdi=zeros(MeanBacLife,Ncells);
KdiFC=zeros(MeanBacLife,Ncells);
KdiR1=zeros(MaxBacLifed,Ncells);
KdiR2=zeros(MaxBacLifed,Ncells);
KdiR3=zeros(MaxBacLifed,Ncells);

Kdx=zeros(MeanBacLife,Ncells);
KdxR1=zeros(MeanBacLifed,Ncells);
KdxR2=zeros(MeanBacLifed,Ncells);
KdxR3=zeros(MeanBacLifed,Ncells);

Kdy=zeros(MeanBacLife,Ncells);
KdyR1=zeros(MeanBacLifed,Ncells);
KdyR2=zeros(MeanBacLifed,Ncells);
KdyR3=zeros(MeanBacLifed,Ncells);

% Use interpolation to make every cell life the same length
% (time synchronisation).

for i=1:Ncells
    %First Channel
    %Spot Integrated Intensities
    S{i}.x(:,1)=imresize(D{i}.x(:,6),[MeanBacLife 1],'bilinear');
     S{i}.xR1(:,1)=imresize(D{i}.xR1(:,6),[MeanBacLife 1],'bilinear');
     S{i}.xR2(:,1)=imresize(D{i}.xR2(:,6),[MeanBacLife 1],'bilinear');
     S{i}.xR3(:,1)=imresize(D{i}.xR3(:,6),[MeanBacLife 1],'bilinear');
    %full cell intensity
    S{i}.x(:,7)=imresize(D{i}.x(:,7),[MeanBacLife 1],'bilinear');
    % x-position
    S{i}.x(:,2)=imresize(D{i}.XNorm(:,2),[MeanBacLife 1],'bilinear');
    S{i}.xR1(:,2)=imresize(D{i}.XNormR1(:,2),[MeanBacLife 1],'bilinear');
    S{i}.xR2(:,2)=imresize(D{i}.XNormR2(:,2),[MeanBacLife 1],'bilinear');
    S{i}.xR3(:,2)=imresize(D{i}.XNormR3(:,2),[MeanBacLife 1],'bilinear');
    % y-position
    S{i}.x(:,4)=imresize(D{i}.XNorm(:,4),[MeanBacLife 1],'bilinear');
    S{i}.xR1(:,4)=imresize(D{i}.XNormR1(:,4),[MeanBacLife 1],'bilinear');
    S{i}.xR2(:,4)=imresize(D{i}.XNormR2(:,4),[MeanBacLife 1],'bilinear');
    S{i}.xR3(:,4)=imresize(D{i}.XNormR3(:,4),[MeanBacLife 1],'bilinear');
    
    % Second Channel
    % spot integrated intensities
    Sd{i}.x(:,1)=imresize(d{i}.x(:,6),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR1(:,1)=imresize(d{i}.xR1(:,6),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR2(:,1)=imresize(d{i}.xR2(:,6),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR3(:,1)=imresize(d{i}.xR3(:,6),[MeanBacLifed 1],'bilinear');
    % full cell integrated intensity
    Sd{i}.x(:,7)=imresize(d{i}.x(:,7),[MeanBacLifed 1],'bilinear');
    % x-position 
    Sd{i}.x(:,2)=imresize(d{i}.XNorm(:,2),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR1(:,2)=imresize(d{i}.XNormR1(:,2),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR2(:,2)=imresize(d{i}.XNormR2(:,2),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR3(:,2)=imresize(d{i}.XNormR3(:,2),[MeanBacLifed 1],'bilinear');
    % y-position
    Sd{i}.x(:,4)=imresize(d{i}.XNorm(:,4),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR1(:,4)=imresize(d{i}.XNormR1(:,4),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR2(:,4)=imresize(d{i}.XNormR2(:,4),[MeanBacLifed 1],'bilinear');
    Sd{i}.xR3(:,4)=imresize(d{i}.XNormR3(:,4),[MeanBacLifed 1],'bilinear');
    
end

%Construct K for taking means of elements at same time point

for j=1:MeanBacLife
    for i=1:Ncells
        Ki(j,i)=S{i}.x(j,1);
        KiFC(j,i)=S{i}.x(j,7);
        KiR1(j,i)=S{i}.xR1(j,1);
        KiR2(j,i)=S{i}.xR2(j,1);
        KiR3(j,i)=S{i}.xR3(j,1);

        Kx(j,i)=S{i}.x(j,2);
        KxR1(j,i)=S{i}.xR1(j,2);
        KxR2(j,i)=S{i}.xR2(j,2);
        KxR3(j,i)=S{i}.xR3(j,2);
        
        Ky(j,i)=S{i}.x(j,4);
        KyR1(j,i)=S{i}.xR1(j,4);
        KyR2(j,i)=S{i}.xR2(j,4);
        KyR3(j,i)=S{i}.xR3(j,4);
    end
end

for j=1:MeanBacLifed
    for i=1:Ncells
        Kdi(j,i)=Sd{i}.x(j,1);
        KdiFC(j,i)=Sd{i}.x(j,7);
        KdiR1(j,i)=Sd{i}.xR1(j,1);
        KdiR2(j,i)=Sd{i}.xR2(j,1);
        KdiR3(j,i)=Sd{i}.xR3(j,1);
        
        Kdx(j,i)=Sd{i}.x(j,2);
        KdxR1(j,i)=Sd{i}.xR1(j,2);
        KdxR2(j,i)=Sd{i}.xR2(j,2);
        KdxR3(j,i)=Sd{i}.xR3(j,2);
        
        Kdy(j,i)=Sd{i}.x(j,4);
        KdyR1(j,i)=Sd{i}.xR1(j,4);
        KdyR2(j,i)=Sd{i}.xR2(j,4);
        KdyR3(j,i)=Sd{i}.xR3(j,4);
        
    end
end

% Mean matrices, which contain means at a time point of relevant parameter
% e.g. M(1,1) contains the mean at the first time point of spot integrated intensity
% e.g. M(:,1) contains the means during cell life of spot integrated intensity
% M(:,1) : spot integrated intensity
% M(:,2) : spot x-position
% M(:,3) : spot x-position std
% M(:,4) : spot y-position
% M(:,5) : spot y-position std
% M(:,6) : spot integrated intensity std
% M(:,7) : full cell integrated intensity
% M(:.8) : full cell integrated intensity std

M=zeros(MeanBacLife,8);
MR1=zeros(MeanBacLife,6);
MR2=zeros(MeanBacLife,6);
MR3=zeros(MeanBacLife,6);

Md=zeros(MeanBacLife,8);
MdR1=zeros(MeanBacLife,6);
MdR2=zeros(MeanBacLife,6);
MdR3=zeros(MeanBacLife,6);

for j=1:MeanBacLife
    M(j,1)=mean(Ki(j,:));
    MR1(j,1)=mean(KiR1(j,:));
    MR2(j,1)=mean(KiR2(j,:));
    MR3(j,1)=mean(KiR3(j,:));
    
    M(j,2)=mean(Kx(j,:));
    MR1(j,2)=mean(KxR1(j,:));
    MR2(j,2)=mean(KxR2(j,:));
    MR3(j,2)=mean(KxR3(j,:));
    
    M(j,3)=std(Kx(j,:));
    MR1(j,3)=std(KxR1(j,:));
    MR2(j,3)=std(KxR2(j,:));
    MR3(j,3)=std(KxR3(j,:));
    
    M(j,4)=mean(Ky(j,:));
    MR1(j,4)=mean(KyR1(j,:));
    MR2(j,4)=mean(KyR2(j,:));
    MR3(j,4)=mean(KyR3(j,:));
    
    M(j,5)=std(Ky(j,:));
    MR1(j,5)=std(KyR1(j,:));
    MR2(j,5)=std(KyR2(j,:));
    MR3(j,5)=std(KyR3(j,:));
    
    M(j,6)=std(Ki(j,:));
    MR1(j,6)=std(KiR1(j,:));
    MR2(j,6)=std(KiR2(j,:));
    MR3(j,6)=std(KiR3(j,:));
    
    M(j,7)=mean(KiFC(j,:));
    M(j,8)=std(KiFC(j,:));
end

for j=1:MeanBacLifed    
    Md(j,1)=mean(Kdi(j,:));
    MdR1(j,1)=mean(KdiR1(j,:));
    MdR2(j,1)=mean(KdiR2(j,:));
    MdR3(j,1)=mean(KdiR3(j,:));
    
    Md(j,2)=mean(Kdx(j,:));
    MdR1(j,2)=mean(KdxR1(j,:));
    MdR2(j,2)=mean(KdxR2(j,:));
    MdR3(j,2)=mean(KdxR3(j,:));
    
    Md(j,3)=std(Kdx(j,:));
    MdR1(j,3)=std(KdxR1(j,:));
    MdR2(j,3)=std(KdxR2(j,:));
    MdR3(j,3)=std(KdxR3(j,:));
    
    Md(j,4)=mean(Kdy(j,:));
    MdR1(j,4)=mean(KdyR1(j,:));
    MdR2(j,4)=mean(KdyR2(j,:));
    MdR3(j,4)=mean(KdyR3(j,:));
    
    Md(j,5)=std(Kdy(j,:));
    MdR1(j,5)=std(KdyR1(j,:));
    MdR2(j,5)=std(KdyR2(j,:));
    MdR3(j,5)=std(KdyR3(j,:));
    
    Md(j,6)=std(Kdi(j,:));
    MdR1(j,6)=std(KdiR1(j,:));
    MdR2(j,6)=std(KdiR2(j,:));
    MdR3(j,6)=std(KdiR3(j,:));
    
    Md(j,7)=mean(KdiFC(j,:));
    Md(j,8)=std(KdiFC(j,:));
end

Mratio=M(:,1)./M(:,7);
Mratiostd=std(Mratio);

Mdratio=Md(:,1)./Md(:,7);
Mdratiostd=std(Mdratio);

xcc=linspace(0,1,MeanBacLife);
xccd=linspace(0,1,MeanBacLifed);

toc

%% Something else

for i=1:Ncells
    Dtotal{i}=D{i}.x(:,6)+D{i}.xR1(:,6);%+D{i}.xR2(:,6)+D{i}.xR3(:,6);
    dtotal{i}=d{i}.x(:,6)+d{i}.xR1(:,6)+d{i}.xR2(:,6);%+d{i}.xR3(:,6);
    
    Dtotsize(i)=size(Dtotal{i},1);
    dtotsize(i)=size(dtotal{i},1);
end

MaxSize=max(Dtotsize);

for i=1:Ncells
    Dtotal{i}(MaxSize+1)=0;
    dtotal{i}(MaxSize+1)=0;
    
    Dtotalz(:,i)=Dtotal{i};
    dtotalz(:,i)=dtotal{i};
end


for j=1:MaxSize+1

nD=Dtotalz(j,:);
nD(nD==0)=[];
LnD=length(nD);
if LnD>0
MeanDtotalz(j,1)=sum(Dtotalz(j,:))./LnD;
else
MeanDtotalz(j,1)=0;
end

nd=dtotalz(j,:);
nd(nd==0)=[];
Lnd=length(nd);
if Lnd>0
Meandtotalz(j,1)=sum(dtotalz(j,:))./Lnd;
else
Meandtotalz(j,1)=0;
end

end

%% Test

hold on
plot(d{1}.x(:,6),'b');
plot(d{1}.xR1(:,6),'r');
hold off
