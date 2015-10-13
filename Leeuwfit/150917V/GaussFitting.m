%% 2D Gaussian fit to spot -- Gaussplosion!!
%
% This code is written to identify multiple spots within the cell. The
% input are the Bacpic images, as a result of the KymoCode, of seperate 
% cells over time. 
%
% written by Roy de Leeuw - 2015
%
% output : x{j}(i,:) = [A,x,sx,y,sy,0]
%
% j : The spot number. Each image could have multiple spots. Stored in cell
% i : The image number. Each cell is image stack as function of time.
% A : Amplitude of the spot.
% x : x-coordinate of the spot. with x = 0 being the cell pole.
% sx : std in x of the spot.
% y : y-coordinate of the spot. 
% sy : std in y of the spot.
%
%% 

clear all
close all
clc

expno='001_DnaN_TUS_dif_30122014_difsignal';
initval=A001_Images_Set_Experiment(expno); %define your paths and files

%% Inputs 

Bac='Fluo0Chan01Bac0003';
BacNum=1;

Mainfolder=strcat(initval.basepath,'StacksLong/dif/');
Stackpth=strcat(Mainfolder,Bac);
d1{1}=readtimeseries(strcat(Stackpth,'/',Bac,'Im'),'tif'); %read zstack
data=dip_array(d1{1}); %turn into uint16 array

%% Ze Defs

Nspots=3; % NUMBER OF SPOTS 

Zsize=size(data,3); XSize=zeros(size(data,2),1);
YSize=zeros(size(data,1),1); %Ysize
ydata=cell(Zsize,1); 
Ydata=cell(Zsize,Nspots); Ydata_X=zeros(Nspots,1); Ydata_Y=zeros(Nspots,1);
ydatacrpd=cell(Zsize,1); ydatacrpdR1=cell(Zsize,1); 
Ampguess=zeros(Zsize,1); Case=cell(Nspots,1); 
x0=cell(Nspots,1); x=cell(Nspots,1); X=cell(Nspots,1); 
Y=cell(Nspots,1); Yg=zeros(Zsize,Nspots); xdata=cell(2,Nspots);
Size=cell(Zsize,Nspots); 

lb = [0,0,0,0,0]; %lower bounds for fit

% fit using levenberg-marquardt algorithm
OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

%% Load Data
% parameters : [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]

% Define some parameters 

SA=3; % pixel step from maximum amplitude 
Sx=3; Sy=3; % guess of sigmaX of spots
Px=3; Py=3; % for optimal fitting add black 'pads' around detected spot
Bs=2; Bss=2*Bs+1; % this defines size of black sqaure
lob=6; upb=14; chanthickness=7; % define the boundaries of the channel

tic 

for i=1:Zsize
    
    %Remove the x,y-Pads
    ydata{i}=double(data(:,:,i));
    [ydatacrpd{i},~]=Crop_Image(ydata{i});
    ydatacrpdR1{i}=ydatacrpd{i};
    
    for j=1:Nspots
    
    % make guesses for spot position
    
    % output = [guessvector,case,newcutoutimage,imageusedforfitting]
    % guessvector = [ampl,x,sigx,y,sigy]
    % case = 2 (spot out of channel), -1 (left edge), 0 (middle), 1 (right
    % edge)
    % imageusedforfitting = only the spot is used for fitting (cut out)
    
    [x0{j}(i,:),Case{j}(i),ydatacrpdR1{i},Ydata{i,j},Size{i,j},Yg(i,j)]= ... 
        ... 
        Spotter(ydatacrpdR1{i},SA,Sx,Sy,Px,Py,Bs,lob,upb);
    
    % preparing framework for gaussian fitting
    [Ydata_X(j),Ydata_Y(j)]=size(Ydata{i,j});
    [X{j},Y{j}] =  meshgrid(linspace(1,Ydata_Y(j),Ydata_Y(j)), ...
        linspace(1,Ydata_X(j),Ydata_X(j)));
    xdata{1,j} = X{j};
    xdata{2,j} = Y{j};
    
    % define upper boundary (Size{i,j}(2) is the X size of the crpd-image
    ub = [realmax('double'),Size{i,j}(2),(Size{i,j}(2))^2,Size{i,j}(2), ...
        (Size{i,j}(2))^2]; %upper bounds
    
    % Do the Gaussian fitting
    
    [x{j}(i,:),resnorm,residual,exitflag] = lsqcurvefit(@GaussPlosFunc, ...
       x0{j}(i,:),xdata(:,j),Ydata{i,j},[],[],OPTIONS);
    
    end
end

%% Translate back to original coords

% use XNorm for normalization to position w.r.t. cell

XNorm=x;

for i=1:Zsize
    for j=1:Nspots
        
    if Case{j}(i)==-1
        x{j}(i,2)=x{j}(i,2)-Px-1;
        x{j}(i,4)=x{j}(i,4)+(Yg(i,j)-SA)-1;
    elseif Case{j}(i)==0
        x{j}(i,2)=x{j}(i,2)-Px+(x0{j}(i,2)-SA)-1;
        x{j}(i,4)=x{j}(i,4)+(Yg(i,j)-SA)-1;
    elseif Case{j}(i)==1
        x{j}(i,2)=x{j}(i,2)-Px-1+(Size{i,j}(2)-2*SA);
        x{j}(i,4)=x{j}(i,4)+(Yg(i,j)-SA)-1;
    elseif Case{j}(i)==2
        x{j}(i,1:5)=0;
        XNorm{j}(i,1:5)=0;
    end

    XNorm{j}(i,2)=x{j}(i,2)/Size{i,j}(2);
    
    XNorm{j}(i,4)=(1-x{j}(i,4)-lob)/chanthickness; 
    % 1- in above, to make sure that 0.9 is upper! part of channel in image
    
    end 
end


%% Integrated Intensity Calculation
II=cell(5,5,Zsize);
FCII=cell(1,Zsize);


for i=1:Zsize
    
    % Full Cell Integrated Intensity
    FCII{i}=ydatacrpd{i}(lob:upb,1:Size{i,j}(2));
    
    for j=1:Nspots
        
    x{j}(i,7)=sum(sum(FCII{i}));
    
        % shift x when it is too close to the left border of image
        if round(x{j}(i,2))==0 || round(x{j}(i,2))==1  
        BorderCase=-1;
        % shift x when it is too close to the right border of image
        elseif round(x{j}(i,2))==Size{i,j}(2) || round(x{j}(i,2))==Size{i,j}(2)-1 
        x{j}(i,2)=Size{i,j}(2)-2;
        BorderCase=1;
        else 
        BorderCase=0;
        end
    
    % SPOT around centroid
    % If sigmaX == 2 or sigmaX == 1
    
    if round(x{j}(i,3))==2 || round(x{j}(i,3))==1 && BorderCase ==0;
    
    for n=[-1 0 1] 
        if round(x{j}(i,5))==2 || round(x{j}(i,5))==1
        for k=[-1 0 1]
    II{n+2,k+2,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+n);
        end
        end
        if round(x(i,5))==3 || round(x(i,5))==4
        for k=[-2 -1 0 1 2]
    II{n+2,k+3,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+n);
        end
        end
    end
    else 
   II{1:3,1:5,i}=0;
    end
    
        % shift x when it is too close to the left border
        if round(x(i,2))==0 || round(x(i,2))==1 || round(x(i,2))==2 
        BorderCase2=-1;
        elseif round(x(i,2))==XSize(i) || round(x(i,2))==XSize(i)-1 || round(x(i,2))==XSize(i)-2 
        % shift x when it is too close to the right border
        BorderCase2=1;
        else 
        BorderCase2=0;
        end
        
    if round(x(i,3))==3 || round(x(i,3))==4 && BorderCase2 == 0;
        

        
        
    for j=[-2 -1 0 1 2] 
        if round(x(i,5))==2 || round(x(i,5))==1
        for k=[-1 0 1]
    II{j+3,k+2,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
        if round(x(i,5))==3 || round(x(i,5))==4
        for k=[-2 -1 0 1 2]
    II{j+3,k+3,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
    end
    end
    
        SII=[II{:,:,i}];
        x(i,6)=sum(SII(:));
    end
end

toc


%% plot

for i=1
% 
figure(1)
hold on
imagesc(ydatacrpd{i})
plot(x(i,2),x(i,4),'r+')
hold off
axis([0 XSize(i)+1 0 size(ydatacrpd{i},1)+1])

xdatafit = linspace(-XSize(i),XSize(i)*2,300);
hdatafit = x(i,1)*exp(-(xdatafit-x(i,2)).^2/(2*x(i,3)^2));
vdatafit = x(i,1)*exp(-(xdatafit-x(i,4)).^2/(2*x(i,5)^2));
% 
figure(2)
hold on
plot(ydatacrpd{i}(round(x(i,4)),:),'ob');
plot(ydatacrpd{i}(:,round(x(i,2))),'or');
plot(xdatafit,hdatafit,'-b',xdatafit,vdatafit,'-r',x(i,2),x(i,6),'g+')
axis([-5 30 0 1000])
hold off
legend('Horizontal','Vertical')

figure(3)
hold on
imagesc(ydatacrpdR1{i})
plot(xR1(i,2),xR1(i,4),'r+')
hold off
axis([0 XSize(i)+1 0 size(ydatacrpdR1{i},1)+1])

hdatafitR1 = xR1(i,1)*exp(-(xdatafit-xR1(i,2)).^2/(2*xR1(i,3)^2));
vdatafitR1 = xR1(i,1)*exp(-(xdatafit-xR1(i,4)).^2/(2*xR1(i,5)^2));
% 
figure(4)
hold on
plot(ydatacrpdR1{i}(round(xR1(i,4)),:),'ob');
plot(ydatacrpdR1{i}(:,round(xR1(i,2))),'or');
plot(xdatafit,hdatafitR1,'-b',xdatafit,vdatafitR1,'-r')
axis([-5 30 0 2000])
hold off
legend('Horizontal','Vertical')
% 
figure(5)
hold on
imagesc(ydatacrpdR2{i})
plot(xR2(i,2),xR2(i,4),'r+')
hold off
axis([0 XSize(i)+1 0 size(ydatacrpdR2{i},1)+1])

hdatafitR2 = xR2(i,1)*exp(-(xdatafit-xR2(i,2)).^2/(2*xR2(i,3)^2));
vdatafitR2 = xR2(i,1)*exp(-(xdatafit-xR2(i,4)).^2/(2*xR2(i,5)^2));

figure(6)
hold on
plot(ydatacrpdR2{i}(round(xR2(i,4)),:),'ob');
plot(ydatacrpdR2{i}(:,round(xR2(i,2))),'or');
plot(xdatafit,hdatafitR2,'-b',xdatafit,vdatafitR2,'-r')
axis([-5 30 0 2000])
hold off
legend('Horizontal','Vertical')

figure(7)
hold on
imagesc(ydatacrpdR3{i})
plot(xR3(i,2),xR3(i,4),'r+')
hold off
axis([0 XSize(i)+1 0 size(ydatacrpdR3{i},1)+1])

hdatafitR3 = xR3(i,1)*exp(-(xdatafit-xR3(i,2)).^2/(2*xR3(i,3)^2));
vdatafitR3 = xR3(i,1)*exp(-(xdatafit-xR3(i,4)).^2/(2*xR3(i,5)^2));

figure(8)
hold on
plot(ydatacrpdR3{i}(round(xR3(i,4)),:),'ob');
plot(ydatacrpdR3{i}(:,round(xR3(i,2))),'or');
plot(xdatafit,hdatafitR3,'-b',xdatafit,vdatafitR3,'-r')
axis([-5 30 0 2000])
hold off
legend('Horizontal','Vertical')

end

%% Save results
 save(strcat(Mainfolder,'DataMULTI/',num2str(BacNum)),'x','xR1','xR2','xR3','XNorm',...
     'XNormR1','XNormR2','XNormR3');