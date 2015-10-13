%% 2D Gaussian fit to spot - RDL 2015
%output : x = [Amp,xcoord,sigmax,ycoord,sigmay,0]
%warning('off')

clear all
close all
clc

initval.basepath='/Users/rleeuw/Work/DataAnalysis/201505_TUSdifDnaN_Montage/Leeuwfit/OriZ_Version/';

%% Inputs 

for z=1
    
BacNum=z;
Bac=num2str(BacNum);

Mainfolder=strcat(initval.basepath,'dif/');
Stackpth=strcat(Mainfolder,Bac,'/');
%d1{1}=readtimeseries(strcat(Stackpth,Bac),'tif'); %read zstack
d1{1}=imread(strcat(Stackpth,Bac));
%data=dip_array(d1{1}); %turn into uint16 array
data=d1{1};

Zsize=size(data,3);
XSize=zeros(size(data,2),1);
YSize=zeros(size(data,1),1);%Ysize
ydata=cell(Zsize,1);%variable that seperate images of stack in cell
Ydata=cell(Zsize,1); YdataR1=cell(Zsize,1); YdataR2=cell(Zsize,1); YdataR3=cell(Zsize,1);
ydatacrpd=cell(Zsize,1); ydatacrpdR1=cell(Zsize,1); ydatacrpdR2=cell(Zsize,1); ydatacrpdR3=cell(Zsize,1);
Xg1=cell(Zsize,1); Xg2=cell(Zsize,1); Xg3=cell(Zsize,1); Xg4=cell(Zsize,1);
Xg5=cell(Zsize,1); Xg6=cell(Zsize,1); Xg7=cell(Zsize,1); Xg8=cell(Zsize,1);
Yg1=cell(Zsize,1); Yg2=cell(Zsize,1); Yg3=cell(Zsize,1); Yg4=cell(Zsize,1);
Yg5=cell(Zsize,1); Yg6=cell(Zsize,1); Yg7=cell(Zsize,1); Yg8=cell(Zsize,1);

Ampguess=zeros(Zsize,1); AmpguessR1=zeros(Zsize,1); AmpguessR2=zeros(Zsize,1);
AmpguessR3=zeros(Zsize,1);
Case=zeros(Zsize,1); CaseR1=zeros(Zsize,1); CaseR2=zeros(Zsize,1); CaseR3=zeros(Zsize,1);

x0=zeros(Zsize,5); x0R1=zeros(Zsize,5); x0R2=zeros(Zsize,5); x0R3=zeros(Zsize,5);
x=zeros(Zsize,5); xR1=zeros(Zsize,5); xR2=zeros(Zsize,5); xR3=zeros(Zsize,5);

XNorm=zeros(Zsize,5); XNormR1=zeros(Zsize,5); XNormR2=zeros(Zsize,5); XNormR3=zeros(Zsize,5);

II=zeros(1,Zsize); IIR1=zeros(1,Zsize); IIR2=zeros(1,Zsize); IIR3=zeros(1,Zsize);
FCII=cell(1,Zsize);

ydatacrpdII=cell(1,Zsize); ydatacrpdIIR1=cell(1,Zsize); ydatacrpdIIR2=cell(1,Zsize); ydatacrpdIIR3=cell(1,Zsize); 

%% Load Data
% parameters : [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
SfA=4; SigX=5; SigY=5; padx=5; pady=5; 
Blackstep=1.5; Blackscratch=2*Blackstep+1;

% Define upper and lower boundary of the channel to filter localizations
upb=18; UB=upb;
lob=6; LB=lob;
chanl=upb-lob;

j=0;
tic

for i=1:Zsize
j=j+1
ydata{i}=double(data(:,:,i));
[ydatacrpd{i},~]=Crop_Image(ydata{i}); %Remove the x,y-Pads of bacpics

% First spot from here -------------------------------------------------
ydatacrpd{i}=[zeros(5,size(ydatacrpd{i},2)) ; ydatacrpd{i} ; zeros(5,size(ydatacrpd{i},2))];
[YSize(i),XSize(i)]=size(ydatacrpd{i});
[Amp,I]=max(ydatacrpd{i});
[Amp2,I2]=max(Amp);
Ampguess(i)=max(Amp2);
Xg1{i}=I2;
Yg1{i}=I(I2);

ydatacrpdR1{i}=ydatacrpd{i};

% Following part is to make an image around the spot of interest 

if (Xg1{i}(1)-SfA)<=0
Ydata{i}=ydatacrpd{i}((Yg1{i}(1)-SfA):(Yg1{i}(1)+SfA),1:2*SfA+1);
Ydata{i}=[zeros(2*SfA+1,padx) Ydata{i} ... 
    zeros(2*SfA+1,padx)];
Case(i)=-1;
% ydatacrpdR1{i}((Yg1{i}(1)-Blackstep):(Yg1{i}(1)+Blackstep),1:Blackscratch)=0;
elseif (Xg1{i}(1)+SfA)>XSize(i)
Ydata{i}=ydatacrpd{i}((Yg1{i}(1)-SfA):(Yg1{i}(1)+SfA),XSize(i)-2*SfA:XSize(i));
Ydata{i}=[zeros(2*SfA+1,padx) Ydata{i} ... 
    zeros(2*SfA+1,padx)];
% ydatacrpdR1{i}((Yg1{i}(1)-Blackstep):(Yg1{i}(1)+Blackstep),XSize(i)-Blackscratch+1:XSize(i))=0;
Case(i)=1;
else 
Ydata{i}=ydatacrpd{i}((Yg1{i}(1)-SfA):(Yg1{i}(1)+SfA),(Xg1{i}(1)-SfA):(Xg1{i}(1)+SfA));
Ydata{i}=[zeros(2*SfA+1,padx) Ydata{i} ... 
    zeros(2*SfA+1,padx)];
% ydatacrpdR1{i}((Yg1{i}(1)-Blackstep):(Yg1{i}(1)+Blackstep),(Xg1{i}(1)-Blackstep):(Xg1{i}(1)+Blackstep))=0;
Case(i)=0;
end

% now we estimate location of spot in new "spot-only" image, and use these
% to guide the fit.

[Dummy,I3]=max(Ydata{i});
[~,I4]=max(Dummy);
Xg2{i}=I4;
Yg2{i}=I3(I4);
[Ydata_Xl,Ydata_Yl]=size(Ydata{i});
x0(i,:)=[Ampguess(i),Xg2{i}(1),SigX,Yg2{i}(1),SigY];% initial value for fit

% FITTING FUNCTION FIRST SPOT ---------------
[X,Y] =  meshgrid(linspace(1,Ydata_Yl,Ydata_Yl),linspace(1,Ydata_Xl,Ydata_Xl));
xdata = cell(2,1);
xdata{1} = X;
xdata{2} = Y;

lb = [0,0,0,0,0];
ub = [realmax('double'),XSize(i),(XSize(i))^2,XSize(i),(XSize(i))^2];

[x(i,:),resnorm,residual,exitflag,~] = lsqcurvefit(@GaussPlosFunc,x0(i,:),xdata,Ydata{i},lb,ub);

% first we translate results to original coordinates:

if Case(i)==-1
        x(i,2)=x(i,2)-padx-1;
elseif Case(i)==0
        x(i,2)=x(i,2)-padx+(Xg1{i}(1)-SfA)-1;
elseif Case(i)==1
        x(i,2)=x(i,2)-padx-1+(XSize(i)-2*SfA);
end

% store normalized data in XNorm, but first define like this.
XNorm(i,:)=x(i,:);
x(i,4)=x(i,4)+(Yg1{i}(1)-SfA)-1;
XNorm(i,2)=x(i,2)/XSize(i);
XNorm(i,4)=(x(i,4)-lob)/chanl;

% now we use the values from fit to determine how much we cut out of the
% original image. Given is x = [amp, x, sigx, y, sigy]. 
siggarX=round(1.5*x(i,3));
siggarY=round(1.5*x(i,5));

RX=round(x(i,2));
RY=round(x(i,4));

if siggarY>3
    siggarY=3;
end

if (RX-siggarX)<=0
ydatacrpdR1{i}((RY-siggarY):(RY+siggarY),1:2*siggarX+1)=0;
elseif (RX+siggarX)>XSize(i);
ydatacrpdR1{i}((RY-siggarY):(RY+siggarY),XSize(i)-2*siggarX+2:XSize(i))=0;   
else 
ydatacrpdR1{i}((RY-siggarY):(RY+siggarY),(RX-siggarX):(RX+siggarX))=0;    
end

% Second Spot From here !!! ---------------------------------------------

[AmpR1,IR1]=max(ydatacrpdR1{i});
[Amp2R1,I2R1]=max(AmpR1);
AmpguessR1(i)=max(Amp2R1);
Xg3{i}=I2R1;
Yg3{i}=IR1(I2R1);

ydatacrpdR2{i}=ydatacrpdR1{i};

if (Yg3{i}(1))<LB || (Yg3{i}(1))>UB % (Yg3(i)-SfA)<=0 || (Yg3(i)+SfA)>YSize(i)
%Make picture dark if spot is not in channel
YdataR1{i}=zeros(2*SfA+1,2*padx+2*SfA+1);
% ydatacrpdR2{i}(:,:)=0;
CaseR1(i)=2; 
elseif (Xg3{i}(1)-SfA)<=0
YdataR1{i}=ydatacrpdR1{i}((Yg3{i}(1)-SfA):(Yg3{i}(1)+SfA),1:2*SfA+1);
YdataR1{i}=[zeros(2*SfA+1,padx) YdataR1{i} ... 
    zeros(2*SfA+1,padx)];
CaseR1(i)=-1;
% ydatacrpdR2{i}((Yg3{i}(1)-Blackstep):(Yg3{i}(1)+Blackstep),1:Blackscratch)=0;
elseif (Xg3{i}(1)+SfA)>XSize(i)
YdataR1{i}=ydatacrpdR1{i}((Yg3{i}(1)-SfA):(Yg3{i}(1)+SfA),XSize(i)-2*SfA:XSize(i));
YdataR1{i}=[zeros(2*SfA+1,padx) YdataR1{i} ... 
    zeros(2*SfA+1,padx)];
% ydatacrpdR2{i}((Yg3{i}(1)-Blackstep):(Yg3{i}(1)+Blackstep),XSize(i)-Blackscratch+1:XSize(i))=0;
CaseR1(i)=1;
else 
YdataR1{i}=ydatacrpdR1{i}((Yg3{i}(1)-SfA):(Yg3{i}(1)+SfA),(Xg3{i}(1)-SfA):(Xg3{i}(1)+SfA));
YdataR1{i}=[zeros(2*SfA+1,padx) YdataR1{i} ... 
    zeros(2*SfA+1,padx)];
% ydatacrpdR2{i}((Yg3{i}(1)-Blackstep):(Yg3{i}(1)+Blackstep),(Xg3{i}(1)-Blackstep):(Xg3{i}(1)+Blackstep))=0;
CaseR1(i)=0;
end

[DummyR1,I3R1]=max(YdataR1{i});
[~,I4R1]=max(DummyR1);
Xg4{i}=I4R1;
Yg4{i}=I3R1(I4R1);
[YdataR1_Xl,YdataR1_Yl]=size(YdataR1{i});
x0R1(i,:)=[AmpguessR1(i),Xg4{i}(1),SigX,Yg4{i}(1),SigY];% initial value 2nd spot

% FITTING FUNCTION SECOND SPOT ---------------
[X2,Y2] =  meshgrid(linspace(1,YdataR1_Yl,YdataR1_Yl),linspace(1,YdataR1_Xl,YdataR1_Xl));
xdata2 = cell(2,1);
xdata2{1} = X2;
xdata2{2} = Y2;

[xR1(i,:),resnorm2,residual2,exitflag2] = lsqcurvefit(@GaussPlosFunc,x0R1(i,:),xdata2,YdataR1{i},lb,ub);

% first we translate results to original coordinates:
    if CaseR1(i)==-1
        xR1(i,2)=xR1(i,2)-padx-1;
    elseif CaseR1(i)==0
        xR1(i,2)=xR1(i,2)-padx+(Xg3{i}(1)-SfA)-1;
    elseif CaseR1(i)==1
        xR1(i,2)=xR1(i,2)-padx-1+(XSize(i)-2*SfA);
    elseif CaseR1(i)==2
        xR1(i,1:5)=1;
    end
    
% store normalized data in XNorm, but first define like this.
XNormR1(i,:)=xR1(i,:);
xR1(i,4)=xR1(i,4)+(Yg1{i}(1)-SfA)-1;
XNormR1(i,2)=xR1(i,2)/XSize(i);
XNormR1(i,4)=(xR1(i,4)-lob)/chanl;

% make sure that spot is within channel
    if xR1(i,4)>upb || xR1(i,4)<lob
        xR1(i,1:5)=0; 
        XNormR1(i,1:5)=0;
    end
    
% now we use the values from fit to determine how much we cut out of the
% original image. Given is x = [amp, x, sigx, y, sigy]. 
siggarXR1=round(1.5*xR1(i,3));
siggarYR1=round(1.5*xR1(i,5));

RXR1=round(xR1(i,2));
RYR1=round(xR1(i,4));

if siggarYR1>3
    siggarYR1=3;
end

if RYR1<4
    RYR1=4;
end

if (RXR1-siggarXR1)<=0
ydatacrpdR2{i}((RYR1-siggarYR1):(RYR1+siggarYR1),1:2*siggarXR1+1)=0;
elseif (RXR1+siggarXR1)>XSize(i);
ydatacrpdR2{i}((RYR1-siggarYR1):(RYR1+siggarYR1),XSize(i)-2*siggarXR1+2:XSize(i))=0;   
else 
ydatacrpdR2{i}((RYR1-siggarYR1):(RYR1+siggarYR1),(RXR1-siggarXR1):(RXR1+siggarXR1))=0;    
end

% make sure that spot is within channel
    if xR1(i,4)>upb || xR1(i,4)<lob
        xR1(i,1:5)=0; 
        XNormR1(i,1:5)=0;
    end
    
%Third Spot From Here! -------------------------------------------------

[AmpR2,IR2]=max(ydatacrpdR2{i});
[Amp2R2,I2R2]=max(AmpR2);
AmpguessR2(i)=max(Amp2R2);
Xg5{i}=I2R2;
Yg5{i}=IR2(I2R2);

ydatacrpdR3{i}=ydatacrpdR2{i};

if (Yg5{i}(1))<LB || (Yg5{i}(1))>UB %(Yg5(i)-SfA)<=0 || (Yg5(i)+SfA)>YSize(i)
% Make picture dark if spot is not in channel
YdataR2{i}=zeros(2*SfA+1+2*pady,2*padx+2*SfA+1);
% ydatacrpdR3{i}(:,:)=0;
CaseR2(i)=2; 
elseif (Xg5{i}(1)-SfA)<=0
YdataR2{i}=ydatacrpdR2{i}((Yg5{i}(1)-SfA):(Yg5{i}(1)+SfA),1:2*SfA+1);
YdataR2{i}=[zeros(2*SfA+1,padx) YdataR2{i} ...
    zeros(2*SfA+1,padx)];
% ydatacrpdR3{i}((Yg5{i}(1)-Blackstep):(Yg5{i}(1)+Blackstep),1:Blackscratch)=0;
CaseR2(i)=-1;
elseif (Xg5{i}(1)+SfA)>XSize(i) 
YdataR2{i}=ydatacrpdR2{i}((Yg5{i}(1)-SfA):(Yg5{i}(1)+SfA),XSize(i)-2*SfA:XSize(i));
YdataR2{i}=[zeros(2*SfA+1,padx) YdataR2{i} ...
    zeros(2*SfA+1,padx)];
% ydatacrpdR3{i}((Yg5{i}(1)-Blackstep):(Yg5{i}(1)+Blackstep),XSize(i)-Blackscratch+1:XSize(i))=0;
CaseR2(i)=1;
else 
YdataR2{i}=ydatacrpdR1{i}((Yg5{i}(1)-SfA):(Yg5{i}(1)+SfA),(Xg5{i}(1)-SfA):(Xg5{i}(1)+SfA));
YdataR2{i}=[zeros(2*SfA+1,padx) YdataR2{i} ...
    zeros(2*SfA+1,padx)];
% ydatacrpdR3{i}((Yg5{i}(1)-Blackstep):(Yg5{i}(1)+Blackstep),(Xg5{i}(1)-Blackstep):(Xg5{i}(1)+Blackstep))=0;
CaseR2(i)=0;
end

[DummyR2,I3R2]=max(YdataR2{i});
[~,I4R2]=max(DummyR2);
Xg6{i}=I4R2;
Yg6{i}=I3R2(I4R2);
[YdataR2_Xl,YdataR2_Yl]=size(YdataR2{i});
x0R2(i,:)=[AmpguessR2(i),Xg6{i}(1),SigX,Yg6{i}(1),SigY]; % initial value 3rd spot

% FITTING FUNCTION THIRD SPOT ---------------
[X3,Y3] =  meshgrid(linspace(1,YdataR2_Yl,YdataR2_Yl),linspace(1,YdataR2_Xl,YdataR2_Xl));
xdata3 = cell(2,1);
xdata3{1} = X3;
xdata3{2} = Y3;


[xR2(i,:),resnorm2,residual2,exitflag2] = lsqcurvefit(@GaussPlosFunc,x0R2(i,:),xdata3,YdataR2{i},lb,ub);

% first we translate results to original coordinates:
    if CaseR2(i)==-1
        xR2(i,2)=xR2(i,2)-padx-1;
    elseif CaseR2(i)==0
        xR2(i,2)=xR2(i,2)-padx+(Xg5{i}(1)-SfA)-1;
    elseif CaseR2(i)==1
        xR2(i,2)=xR2(i,2)-padx-1+(XSize(i)-2*SfA);
    elseif CaseR2(i)==2
        xR2(i,1:5)=1; 
    end

% store normalized data in XNorm, but first define like this.
XNormR2(i,:)=xR2(i,:);
xR2(i,4)=xR2(i,4)+(Yg1{i}(1)-SfA)-1;
XNormR2(i,2)=xR2(i,2)/XSize(i);
XNormR2(i,4)=(xR2(i,4)-lob)/chanl;
    
% now we use the values from fit to determine how much we cut out of the
% original image. Given is x = [amp, x, sigx, y, sigy]. 
siggarXR2=round(1.5*xR2(i,3));
siggarYR2=round(1.5*xR2(i,5));

RXR2=round(xR2(i,2));
RYR2=round(xR2(i,4));

if siggarYR2>3
    siggarYR2=3;
end

if RYR2<5
    RYR2=5;
end

if (RXR2-siggarXR2)<=0
ydatacrpdR3{i}((RYR2-siggarYR2):(RYR2+siggarYR2),1:2*siggarXR2+1)=0;
elseif (RXR2+siggarXR2)>XSize(i);
ydatacrpdR3{i}((RYR2-siggarYR2):(RYR2+siggarYR2),XSize(i)-2*siggarXR2+2:XSize(i))=0;   
else 
ydatacrpdR3{i}((RYR2-siggarYR2):(RYR2+siggarYR2),(RXR2-siggarXR2):(RXR2+siggarXR2))=0;    
end

% make sure that spot is within channel
    if xR2(i,4)>upb || xR2(i,4)<lob
        xR2(i,1:5)=0; 
        XNormR2(i,1:5)=0;
    end
    
%Fourth Spot from Here!!!  ----------------------------------------------

[AmpR3,IR3]=max(ydatacrpdR3{i});
[Amp2R3,I2R3]=max(AmpR3);
AmpguessR3(i)=max(Amp2R3);
Xg7{i}=I2R3;
Yg7{i}=IR3(I2R3);

if (Yg7{i}(1))<LB || (Yg7{i}(1))>UB
YdataR3{i}=zeros(2*SfA+1+2*pady,2*padx+2*SfA+1);
CaseR3(i)=2;
elseif (Xg7{i}(1)-SfA)<=0 
YdataR3{i}=ydatacrpdR2{i}((Yg7{i}(1)-SfA):(Yg7{i}(1)+SfA),1:2*SfA+1);
YdataR3{i}=[zeros(2*SfA+1,padx) YdataR3{i} ...
    zeros(2*SfA+1,padx)];
CaseR3(i)=-1;
elseif (Xg7{i}(1)+SfA)>XSize(i) 
YdataR3{i}=ydatacrpdR2{i}((Yg7{i}(1)-SfA):(Yg7{i}(1)+SfA),(XSize(i)-2*SfA):XSize(i));
YdataR3{i}=[zeros(2*SfA+1,padx) YdataR3{i} ...
    zeros(2*SfA+1,padx)];
CaseR3(i)=1;
else 
YdataR3{i}=ydatacrpdR1{i}((Yg7{i}(1)-SfA):(Yg7{i}(1)+SfA),(Xg7{i}(1)-SfA):(Xg7{i}(1)+SfA));
YdataR3{i}=[zeros(2*SfA+1,padx) YdataR3{i} ...
    zeros(2*SfA+1,padx)];
CaseR3(i)=0;
end

[DummyR3,I3R3]=max(YdataR3{i});
[~,I4R3]=max(DummyR3);
Xg8{i}=I4R3;
Yg8{i}=I3R3(I4R3);
[YdataR3_Xl,YdataR3_Yl]=size(YdataR3{i});
x0R3(i,:)=[AmpguessR3(i),Xg8{i}(1),SigX,Yg8{i}(1),SigY];

% FITTING FUNCTION FOURTH SPOT -------------------------------------------
[X4,Y4] =  meshgrid(linspace(1,YdataR3_Yl,YdataR3_Yl),linspace(1,YdataR3_Xl,YdataR3_Xl));
xdata4 = cell(2,1);
xdata4{1} = X4;
xdata4{2} = Y4;

[xR3(i,:),resnorm4,residual4,exitflag4] = lsqcurvefit(@GaussPlosFunc,x0R3(i,:),xdata4,YdataR3{i},lb,ub);

% first we translate results to original coordinates:
    if CaseR3(i)==-1
        xR3(i,2)=xR3(i,2)-padx-1;
    elseif CaseR3(i)==0
        xR3(i,2)=xR3(i,2)-padx+(Xg7{i}(1)-SfA)-1;
    elseif CaseR3(i)==1
        xR3(i,2)=xR3(i,2)-padx-1+(XSize(i)-2*SfA);
    elseif CaseR3(i)==2
        xR3(i,1:5)=0;
    end

% store normalized data in XNorm, but first define like this.
XNormR3(i,:)=xR3(i,:);
xR3(i,4)=xR3(i,4)+(Yg1{i}(1)-SfA)-1;
XNormR3(i,2)=xR3(i,2)/XSize(i);
XNormR3(i,4)=(xR3(i,4)-lob)/chanl;

% make sure that spot is within channel
    if xR3(i,4)>upb || xR3(i,4)<lob
        xR3(i,1:5)=0;
        XNormR3(i,1:5)=0;
    end
end

PX=round(x(:,2)); PXR1=round(xR1(:,2)); PXR2=round(xR2(:,2)); PXR3=round(xR3(:,2));
PY=round(x(:,4)); PYR1=round(xR1(:,4)); PYR2=round(xR2(:,4)); PYR3=round(xR3(:,4));

SX=round(1.5*x(:,3)); SXR1=round(1.5*xR1(:,3)); SXR2=round(1.5*xR2(:,3)); SXR3=round(1.5*xR3(:,3));
SY=round(1.5*x(:,5)); SYR1=round(1.5*xR1(:,5)); SYR2=round(1.5*xR2(:,5)); SYR3=round(1.5*xR3(:,5));

% ------------------ INTEGRATED INTENSITY CALCULATIONS----------------- %
% first add pad in case SD of spot is out of the picture for II calcs.
% then shift the P's accordingly
% then sum

    PX=PX+padx; PXR1=PXR1+padx; PXR2=PXR2+padx; PXR3=PXR3+padx; 
    PY=PY+pady; PYR1=PYR1+pady; PYR2=PYR2+pady; PYR3=PYR3+pady; 
    SY(SY>6)=6; SYR1(SYR1>6)=6; SYR2(SYR2>6)=6; SYR3(SYR3>6)=6;
    
for i=1:Zsize

    ydatacrpdII{i}=[zeros(pady,size(ydatacrpd{i},2)+2*padx); zeros(size(ydatacrpd{i},1),padx) ydatacrpd{i} ...
        zeros(size(ydatacrpd{i},1),padx) ;zeros(pady,size(ydatacrpd{i},2)+2*padx)];
    ydatacrpdIIR1{i}=[zeros(pady,size(ydatacrpdR1{i},2)+2*padx); zeros(size(ydatacrpdR1{i},1),padx) ydatacrpdR1{i} ...
        zeros(size(ydatacrpdR1{i},1),padx) ;zeros(pady,size(ydatacrpdR1{i},2)+2*padx)];
    ydatacrpdIIR2{i}=[zeros(pady,size(ydatacrpdR2{i},2)+2*padx); zeros(size(ydatacrpdR2{i},1),padx) ydatacrpdR2{i} ...
        zeros(size(ydatacrpdR2{i},1),padx) ;zeros(pady,size(ydatacrpdR2{i},2)+2*padx)];
    ydatacrpdIIR3{i}=[zeros(pady,size(ydatacrpdR3{i},2)+2*padx); zeros(size(ydatacrpdR3{i},1),padx) ydatacrpdR3{i} ... 
        zeros(size(ydatacrpdR3{i},1),padx) ;zeros(pady,size(ydatacrpdR3{i},2)+2*padx)];


    % Integrated Intensities
    x(i,6)=sum(sum(ydatacrpdII{i}(PY(i)-SY(i):PY(i)+SY(i),PX(i)-SX(i):PX(i)+SX(i))));
    xR1(i,6)=sum(sum(ydatacrpdIIR1{i}(PYR1(i)-SYR1(i):PYR1(i)+SYR1(i),PXR1(i)-SXR1(i):PXR1(i)+SXR1(i))));
    xR2(i,6)=sum(sum(ydatacrpdIIR2{i}(PYR2(i)-SYR2(i):PYR2(i)+SYR2(i),PXR2(i)-SXR2(i):PXR2(i)+SXR2(i))));
    xR3(i,6)=sum(sum(ydatacrpdIIR3{i}(PYR3(i)-SYR3(i):PYR3(i)+SYR3(i),PXR3(i)-SXR3(i):PXR3(i)+SXR3(i))));

    % Full Cell Intensities
    x(i,7)=sum(sum(ydatacrpd{i}(lob:upb,1:XSize(i))));
    
end

toc




%% plot

% for i=15
% 
% figure(1)
% hold on
% imagesc(ydatacrpd{i})
% plot(x(i,2),x(i,4),'r+')
% hold off
% axis([0 XSize(i)+1 0 size(ydatacrpd{i},1)+1])
% 
% xdatafit = linspace(-XSize(i),XSize(i)*2,300);
% hdatafit = x(i,1)*exp(-(xdatafit-x(i,2)).^2/(2*x(i,3)^2));
% vdatafit = x(i,1)*exp(-(xdatafit-x(i,4)).^2/(2*x(i,5)^2));
% 
% figure(2)
% hold on
% plot(ydatacrpd{i}(round(x(i,4)),:),'ob');
% plot(ydatacrpd{i}(:,round(x(i,2))),'or');
% plot(xdatafit,hdatafit,'-b',xdatafit,vdatafit,'-r',x(i,2),x(i,6),'g+')
% axis([0 35 0 9000])
% hold off
% legend('Horizontal','Vertical')
% xlabel('Position (-)')
% ylabel('Intensity Counts (-)')
% set(gca,'FontSize',15)
% 
% 
% figure(3)
% hold on
% imagesc(ydatacrpdR1{i})
% plot(xR1(i,2),xR1(i,4),'r+')
% hold off
% axis([0 XSize(i)+1 0 size(ydatacrpdR1{i},1)+1])
% 
% hdatafitR1 = xR1(i,1)*exp(-(xdatafit-xR1(i,2)).^2/(2*xR1(i,3)^2));
% vdatafitR1 = xR1(i,1)*exp(-(xdatafit-xR1(i,4)).^2/(2*xR1(i,5)^2));
% 
% figure(4)
% hold on
% plot(ydatacrpdR1{i}(round(xR1(i,4)),:),'ob');
% plot(ydatacrpdR1{i}(:,round(xR1(i,2))),'or');
% plot(xdatafit,hdatafitR1,'-b',xdatafit,vdatafitR1,'-r')
% axis([-5 60 0 5000])
% hold off
% legend('Horizontal','Vertical')
% xlabel('Position (-)')
% ylabel('Intensity Counts (-)')
% set(gca,'FontSize',15)
% 
% figure(5)
% hold on
% imagesc(ydatacrpdR2{i})
% plot(xR2(i,2),xR2(i,4),'r+')
% hold off
% axis([0 XSize(i)+1 0 size(ydatacrpdR2{i},1)+1])
% 
% hdatafitR2 = xR2(i,1)*exp(-(xdatafit-xR2(i,2)).^2/(2*xR2(i,3)^2));
% vdatafitR2 = xR2(i,1)*exp(-(xdatafit-xR2(i,4)).^2/(2*xR2(i,5)^2));
% 
% figure(6)
% hold on
% plot(ydatacrpdR2{i}(round(xR2(i,4)),:),'ob');
% plot(ydatacrpdR2{i}(:,round(xR2(i,2))),'or');
% plot(xdatafit,hdatafitR2,'-b',xdatafit,vdatafitR2,'-r')
% axis([0 30 0 2500])
% hold off
% legend('Horizontal','Vertical')
% xlabel('Position (-)')
% ylabel('Intensity Counts (-)')
% set(gca,'FontSize',15)
% 
% figure(7)
% hold on
% imagesc(ydatacrpdR3{i})
% plot(xR3(i,2),xR3(i,4),'r+')
% hold off
% axis([0 XSize(i)+1 0 size(ydatacrpdR3{i},1)+1])
% 
% hdatafitR3 = xR3(i,1)*exp(-(xdatafit-xR3(i,2)).^2/(2*xR3(i,3)^2));
% vdatafitR3 = xR3(i,1)*exp(-(xdatafit-xR3(i,4)).^2/(2*xR3(i,5)^2));
% 
% figure(8)
% hold on
% plot(ydatacrpdR3{i}(round(xR3(i,4)),:),'ob');
% plot(ydatacrpdR3{i}(:,round(xR3(i,2))),'or');
% plot(xdatafit,hdatafitR3,'-b',xdatafit,vdatafitR3,'-r')
% axis([-5 60 0 5000])
% hold off
% legend('Horizontal','Vertical')
% xlabel('Position (-)')
% ylabel('Intensity Counts (-)')
% set(gca,'FontSize',15)
% 
% end
end
%% Save results
 save(strcat(Mainfolder,'DataMULTI/',num2str(BacNum)),'x','xR1','xR2','xR3','XNorm',...
     'XNormR1','XNormR2','XNormR3','XSize');
 disp(strcat(Bac,' done')); 
