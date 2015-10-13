%% 2D Gaussian fit to spot -- Gaussplosion!!
% output : x = [Amp,xcoord,sigmax,ycoord,sigmay,0]
%warning('off')

clear all
close all
clc

expno='001_DnaN_TUS_dif_30122014_TUSsignal';
initval=A001_Images_Set_Experiment(expno); %define your paths and files

%% Inputs 

Bac='Fluo0Chan02Bac0049';
BacNum=18;

Mainfolder=strcat(initval.basepath,'StacksLong/DnaN/');
Stackpth=strcat(Mainfolder,Bac);
d1{1}=readtimeseries(strcat(Stackpth,'/',Bac,'Im'),'tif'); %read zstack
data=dip_array(d1{1}); %turn into uint16 array

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

%% Load Data
% parameters : [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
SfA=3; SigX=3; SigY=3; padx=3; pady=3; 
Blackstep=2; Blackscratch=2*Blackstep+1;

tic
for i=1:Zsize
    
ydata{i}=double(data(:,:,i));
[ydatacrpd{i},~]=Crop_Image(ydata{i}); %Remove the x,y-Pads

% First spot from here -------------------------------------------------

[YSize(i),XSize(i)]=size(ydatacrpd{i});
[Amp,I]=max(ydatacrpd{i});
[Amp2,I2]=max(Amp);
Ampguess(i)=max(Amp2);
Xg1{i}=I2;
Yg1{i}=I(I2);

ydatacrpdR1{i}=ydatacrpd{i};

if (Xg1{i}(1)-SfA)<=0
Ydata{i}=ydatacrpd{i}((Yg1{i}(1)-SfA):(Yg1{i}(1)+SfA),1:2*SfA+1);
Ydata{i}=[zeros(2*SfA+1,padx) Ydata{i} zeros(2*SfA+1,padx)];
Case(i)=-1;
ydatacrpdR1{i}((Yg1{i}(1)-Blackstep):(Yg1{i}(1)+Blackstep),1:Blackscratch)=0;
elseif (Xg1{i}(1)+SfA)>XSize(i)
Ydata{i}=ydatacrpd{i}((Yg1{i}(1)-SfA):(Yg1{i}(1)+SfA),XSize(i)-2*SfA:XSize(i));
Ydata{i}=[zeros(2*SfA+1,padx) Ydata{i} zeros(2*SfA+1,padx)];
ydatacrpdR1{i}((Yg1{i}(1)-Blackstep):(Yg1{i}(1)+Blackstep),XSize(i)-Blackscratch+1:XSize(i))=0;
Case(i)=1;
else 
Ydata{i}=ydatacrpd{i}((Yg1{i}(1)-SfA):(Yg1{i}(1)+SfA),(Xg1{i}(1)-SfA):(Xg1{i}(1)+SfA));
Ydata{i}=[zeros(2*SfA+1,padx) Ydata{i} zeros(2*SfA+1,padx)];
ydatacrpdR1{i}((Yg1{i}(1)-Blackstep):(Yg1{i}(1)+Blackstep),(Xg1{i}(1)-Blackstep):(Xg1{i}(1)+Blackstep))=0;
Case(i)=0;
end

[Dummy,I3]=max(Ydata{i});
[~,I4]=max(Dummy);
Xg2{i}=I4;
Yg2{i}=I3(I4);
[Ydata_Xl,Ydata_Yl]=size(Ydata{i});
x0(i,:)=[Ampguess(i),Xg2{i}(1),SigX,Yg2{i}(1),SigY];% initial value

% Second Spot From here !!! ---------------------------------------------

[AmpR1,IR1]=max(ydatacrpdR1{i});
[Amp2R1,I2R1]=max(AmpR1);
AmpguessR1(i)=max(Amp2R1);
Xg3{i}=I2R1;
Yg3{i}=IR1(I2R1);

ydatacrpdR2{i}=ydatacrpdR1{i};

if (Yg3{i}(1))<6 || (Yg3{i}(1))>14 % (Yg3(i)-SfA)<=0 || (Yg3(i)+SfA)>YSize(i)
%Make picture dark if spot is not in channel
YdataR1{i}=zeros(2*SfA+1,2*padx+2*SfA+1);
ydatacrpdR2{i}(:,:)=0;
CaseR1(i)=2; 
elseif (Xg3{i}(1)-SfA)<=0
YdataR1{i}=ydatacrpdR1{i}((Yg3{i}(1)-SfA):(Yg3{i}(1)+SfA),1:2*SfA+1);
YdataR1{i}=[zeros(2*SfA+1,padx) YdataR1{i} zeros(2*SfA+1,padx)];
CaseR1(i)=-1;
ydatacrpdR2{i}((Yg3{i}(1)-Blackstep):(Yg3{i}(1)+Blackstep),1:Blackscratch)=0;
elseif (Xg3{i}(1)+SfA)>XSize(i)
YdataR1{i}=ydatacrpdR1{i}((Yg3{i}(1)-SfA):(Yg3{i}(1)+SfA),XSize(i)-2*SfA:XSize(i));
YdataR1{i}=[zeros(2*SfA+1,padx) YdataR1{i} zeros(2*SfA+1,padx)];
ydatacrpdR2{i}((Yg3{i}(1)-Blackstep):(Yg3{i}(1)+Blackstep),XSize(i)-Blackscratch+1:XSize(i))=0;
CaseR1(i)=1;
else 
YdataR1{i}=ydatacrpdR1{i}((Yg3{i}(1)-SfA):(Yg3{i}(1)+SfA),(Xg3{i}(1)-SfA):(Xg3{i}(1)+SfA));
YdataR1{i}=[zeros(2*SfA+1,padx) YdataR1{i} zeros(2*SfA+1,padx)];
ydatacrpdR2{i}((Yg3{i}(1)-Blackstep):(Yg3{i}(1)+Blackstep),(Xg3{i}(1)-Blackstep):(Xg3{i}(1)+Blackstep))=0;
CaseR1(i)=0;
end

[DummyR1,I3R1]=max(YdataR1{i});
[~,I4R1]=max(DummyR1);
Xg4{i}=I4R1;
Yg4{i}=I3R1(I4R1);
[YdataR1_Xl,YdataR1_Yl]=size(YdataR1{i});
x0R1(i,:)=[AmpguessR1(i),Xg4{i}(1),SigX,Yg4{i}(1),SigY];% initial value 2nd spot

%Third Spot From Here! -------------------------------------------------

[AmpR2,IR2]=max(ydatacrpdR2{i});
[Amp2R2,I2R2]=max(AmpR2);
AmpguessR2(i)=max(Amp2R2);
Xg5{i}=I2R2;
Yg5{i}=IR2(I2R2);

ydatacrpdR3{i}=ydatacrpdR2{i};

if (Yg5{i}(1))<6 || (Yg5{i}(1))>14 %(Yg5(i)-SfA)<=0 || (Yg5(i)+SfA)>YSize(i)
% Make picture dark if spot is not in channel
YdataR2{i}=zeros(2*SfA+1+2*pady,2*padx+2*SfA+1);
ydatacrpdR3{i}(:,:)=0;
CaseR2(i)=2; 
elseif (Xg5{i}(1)-SfA)<=0
YdataR2{i}=ydatacrpdR2{i}((Yg5{i}(1)-SfA):(Yg5{i}(1)+SfA),1:2*SfA+1);
YdataR2{i}=[zeros(pady,2*padx+size(YdataR2{i},2)); zeros(2*SfA+1,padx) YdataR2{i} ...
    zeros(2*SfA+1,padx); zeros(pady,2*padx+size(YdataR2{i},2))];
ydatacrpdR3{i}((Yg5{i}(1)-Blackstep):(Yg5{i}(1)+Blackstep),1:Blackscratch)=0;
CaseR2(i)=-1;
elseif (Xg5{i}(1)+SfA)>XSize(i) 
YdataR2{i}=ydatacrpdR2{i}((Yg5{i}(1)-SfA):(Yg5{i}(1)+SfA),XSize(i)-2*SfA:XSize(i));
YdataR2{i}=[zeros(pady,2*padx+size(YdataR2{i},2)); zeros(2*SfA+1,padx) YdataR2{i} ...
    zeros(2*SfA+1,padx); zeros(pady,2*padx+size(YdataR2{i},2))];
ydatacrpdR3{i}((Yg5{i}(1)-Blackstep):(Yg5{i}(1)+Blackstep),XSize(i)-Blackscratch+1:XSize(i))=0;
CaseR2(i)=1;
else 
YdataR2{i}=ydatacrpdR1{i}((Yg5{i}(1)-SfA):(Yg5{i}(1)+SfA),(Xg5{i}(1)-SfA):(Xg5{i}(1)+SfA));
YdataR2{i}=[zeros(pady,2*padx+size(YdataR2{i},2)); zeros(2*SfA+1,padx) YdataR2{i} ...
    zeros(2*SfA+1,padx); zeros(pady,2*padx+size(YdataR2{i},2))];
ydatacrpdR3{i}((Yg5{i}(1)-Blackstep):(Yg5{i}(1)+Blackstep),(Xg5{i}(1)-Blackstep):(Xg5{i}(1)+Blackstep))=0;
CaseR2(i)=0;
end

[DummyR2,I3R2]=max(YdataR2{i});
[~,I4R2]=max(DummyR2);
Xg6{i}=I4R2;
Yg6{i}=I3R2(I4R2);
[YdataR2_Xl,YdataR2_Yl]=size(YdataR2{i});
x0R2(i,:)=[AmpguessR2(i),Xg6{i}(1),SigX,Yg6{i}(1),SigY]; % initial value 3rd spot

%Fourth Spot from Here!!!  ----------------------------------------------

[AmpR3,IR3]=max(ydatacrpdR3{i});
[Amp2R3,I2R3]=max(AmpR3);
AmpguessR3(i)=max(Amp2R3);
Xg7{i}=I2R3;
Yg7{i}=IR3(I2R3);

if (Yg7{i}(1))<6 || (Yg7{i}(1))>14
YdataR3{i}=zeros(2*SfA+1+2*pady,2*padx+2*SfA+1);
CaseR3(i)=2;
elseif (Xg7{i}(1)-SfA)<=0 
YdataR3{i}=ydatacrpdR2{i}((Yg7{i}(1)-SfA):(Yg7{i}(1)+SfA),1:2*SfA+1);
YdataR3{i}=[zeros(pady,2*padx+size(YdataR3{i},2)); zeros(2*SfA+1,padx) YdataR3{i} ...
    zeros(2*SfA+1,padx); zeros(pady,2*padx+size(YdataR3{i},2))];
CaseR3(i)=-1;
elseif (Xg7{i}(1)+SfA)>XSize(i) 
YdataR3{i}=ydatacrpdR2{i}((Yg7{i}(1)-SfA):(Yg7{i}(1)+SfA),(XSize(i)-2*SfA):XSize(i));
YdataR3{i}=[zeros(pady,2*padx+size(YdataR3{i},2)); zeros(2*SfA+1,padx) YdataR3{i} ...
    zeros(2*SfA+1,padx); zeros(pady,2*padx+size(YdataR3{i},2))];
CaseR3(i)=1;
else 
YdataR3{i}=ydatacrpdR1{i}((Yg7{i}(1)-SfA):(Yg7{i}(1)+SfA),(Xg7{i}(1)-SfA):(Xg7{i}(1)+SfA));
YdataR3{i}=[zeros(pady,2*padx+size(YdataR3{i},2)); zeros(2*SfA+1,padx) YdataR3{i} ...
    zeros(2*SfA+1,padx); zeros(pady,2*padx+size(YdataR3{i},2))];
CaseR3(i)=0;
end

[DummyR3,I3R3]=max(YdataR3{i});
[~,I4R3]=max(DummyR3);
Xg8{i}=I4R3;
Yg8{i}=I3R3(I4R3);
[YdataR3_Xl,YdataR3_Yl]=size(YdataR3{i});
x0R3(i,:)=[AmpguessR3(i),Xg8{i}(1),SigX,Yg8{i}(1),SigY];

end
%% Inputs
% grid for x data;
[X,Y] =  meshgrid(linspace(1,Ydata_Yl,Ydata_Yl),linspace(1,Ydata_Xl,Ydata_Xl));
xdata = cell(2,1);
xdata{1} = X;
xdata{2} = Y;

[X2,Y2] =  meshgrid(linspace(1,YdataR2_Yl,YdataR2_Yl),linspace(1,YdataR2_Xl,YdataR2_Xl));
xdata2 = cell(2,1);
xdata2{1} = X2;
xdata2{2} = Y2;

[X3,Y3] =  meshgrid(linspace(1,YdataR3_Yl,YdataR3_Yl),linspace(1,YdataR3_Xl,YdataR3_Xl));
xdata3 = cell(2,1);
xdata3{1} = X3;
xdata3{2} = Y3;

%% Fit

lb = [0,0,0,0,0]; %lower bounds
for i=1:Zsize
ub = [realmax('double'),XSize(i),(XSize(i))^2,XSize(i),(XSize(i))^2]; %upper bounds
[x(i,:),resnorm,residual,exitflag] = lsqcurvefit(@GaussPlosFunc,x0(i,:),xdata,Ydata{i},lb,ub);
[xR1(i,:),resnorm2,residual2,exitflag2] = lsqcurvefit(@GaussPlosFunc,x0R1(i,:),xdata,YdataR1{i},lb,ub);
[xR2(i,:),resnorm3,residual3,exitflag3] = lsqcurvefit(@GaussPlosFunc,x0R2(i,:),xdata2,YdataR2{i},lb,ub);
[xR3(i,:),resnorm4,residual4,exitflag4] = lsqcurvefit(@GaussPlosFunc,x0R3(i,:),xdata3,YdataR3{i},lb,ub);
end

%% Translate back to original coords

% use XNorm for normalization to position w.r.t. cell

XNorm=x;
XNormR1=xR1;
XNormR2=xR2;
XNormR3=xR3;

% Define upper and lower boundary of the channel to filter localizations
upb=14;
lob=7;
chanl=7;

for i=1:Zsize
    if Case(i)==-1
        x(i,2)=x(i,2)-padx-1;
    elseif Case(i)==0
        x(i,2)=x(i,2)-padx+(Xg1{i}(1)-SfA)-1;
    elseif Case(i)==1
        x(i,2)=x(i,2)-padx-1+(XSize(i)-2*SfA);
    end
    
    x(i,4)=x(i,4)+(Yg1{i}(1)-SfA)-1;
    XNorm(i,2)=x(i,2)/XSize(i);
    XNorm(i,4)=(x(i,4)-lob)/chanl;
    
    if CaseR1(i)==-1
        xR1(i,2)=xR1(i,2)-padx-1;
    elseif CaseR1(i)==0
        xR1(i,2)=xR1(i,2)-padx+(Xg3{i}(1)-SfA)-1;
    elseif CaseR1(i)==1
        xR1(i,2)=xR1(i,2)-padx-1+(XSize(i)-2*SfA);
    elseif CaseR1(i)==2
        xR1(i,1:5)=1;
    end
    if CaseR2(i)==-1
        xR2(i,2)=xR2(i,2)-padx-1;
    elseif CaseR2(i)==0
        xR2(i,2)=xR2(i,2)-padx+(Xg5{i}(1)-SfA)-1;
    elseif CaseR2(i)==1
        xR2(i,2)=xR2(i,2)-padx-1+(XSize(i)-2*SfA);
    elseif CaseR2(i)==2
        xR2(i,1:5)=1; 
    end
    if CaseR3(i)==-1
        xR3(i,2)=xR3(i,2)-padx-1;
    elseif CaseR3(i)==0
        xR3(i,2)=xR3(i,2)-padx+(Xg7{i}(1)-SfA)-1;
    elseif CaseR3(i)==1
        xR3(i,2)=xR3(i,2)-padx-1+(XSize(i)-2*SfA);
    elseif CaseR3(i)==2
        xR3(i,1:5)=0;
    end



xR1(i,4)=xR1(i,4)+(Yg3{i}(1)-SfA)-1;
XNormR1(i,4)=(xR1(i,4)-lob)/chanl;
XNormR1(i,2)=xR1(i,2)/XSize(i);

xR2(i,4)=xR2(i,4)+(Yg5{i}(1)-SfA)-1-pady;
XNormR2(i,4)=(xR2(i,4)-lob)/chanl;
XNormR2(i,2)=xR2(i,2)/XSize(i);
        
xR3(i,4)=xR3(i,4)+(Yg7{i}(1)-SfA)-1-pady;
XNormR3(i,2)=xR3(i,2)/XSize(i);
XNormR3(i,4)=(xR3(i,4)-lob)/chanl;  

    if xR1(i,4)>upb || xR1(i,4)<lob
        xR1(i,1:5)=0; 
        XNormR1(i,1:5)=0;
    end
    if xR2(i,4)>upb || xR2(i,4)<lob
        xR2(i,1:5)=0; 
        XNormR2(i,1:5)=0;
    end
    if xR3(i,4)>upb || xR3(i,4)<lob
        xR3(i,1:5)=0;
        XNormR3(i,1:5)=0;
    end
end
    
%% Integrated Intensity Calculation
II=cell(5,5,Zsize);
II1=cell(5,5,Zsize);
II2=cell(5,5,Zsize);
II3=cell(5,5,Zsize);
FCII=cell(1,Zsize);


for i=1:Zsize
    % Full Cell Integrated Intensity
    FCII{i}=ydatacrpd{i}(6:14,1:XSize(i));
    x(i,7)=sum(sum(FCII{i}));
    
    % SPOT around centroid
    
    if round(x(i,3))==2 || round(x(i,3))==1
        
        % shift x when it is too close to the left border
        if round(x(i,2))==0 || round(x(i,2))==1  
        x(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if round(x(i,2))==XSize(i) || round(x(i,2))==XSize(i)-1 
        x(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(x(i,5))==2 || round(x(i,5))==1
        for k=[-1 0 1]
    II{j+2,k+2,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
        if round(x(i,5))==3 || round(x(i,5))==4
        for k=[-2 -1 0 1 2]
    II{j+2,k+3,i}=ydatacrpd{i}(round(x(i,4))+k,round(x(i,2))+j);
        end
        end
    end
    end
    
    if round(x(i,3))==3 || round(x(i,3))==4
        
        % shift x when it is too close to the left border
        if round(x(i,2))==0 || round(x(i,2))==1 || round(x(i,2))==2 
        x(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(x(i,2))==XSize(i) || round(x(i,2))==XSize(i)-1 || round(x(i,2))==XSize(i)-2 
        x(i,2)=XSize(i)-3;
        end
        
        
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

%---------- Second Spot ---------------------------------------------------
    if round(xR1(i,3))==2 || round(xR1(i,3))==1
        
        % shift x when it is too close to the left border
        if round(xR1(i,2))==0 || round(xR1(i,2))==1  
        xR1(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if round(xR1(i,2))==XSize(i) || round(xR1(i,2))==XSize(i)-1 
        xR1(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(xR1(i,5))==2 || round(xR1(i,5))==1
        for k=[-1 0 1]
    II1{j+2,k+2,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
        if round(xR1(i,5))==3 || round(xR1(i,5))==4
        for k=[-2 -1 0 1 2]
    II1{j+2,k+3,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
    end
    end
    
    if round(xR1(i,3))==3 || round(xR1(i,3))==4
        
        % shift x when it is too close to the left border
        if round(xR1(i,2))==0 || round(xR1(i,2))==1 || round(xR1(i,2))==2 
        xR1(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(xR1(i,2))==XSize(i) || round(xR1(i,2))==XSize(i)-1 || round(xR1(i,2))==XSize(i)-2 
        xR1(i,2)=XSize(i)-3;
        end
        
        
    for j=[-2 -1 0 1 2] 
        if round(xR1(i,5))==2 || round(xR1(i,5))==1
        for k=[-1 0 1]
    II1{j+3,k+2,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
        if round(xR1(i,5))==3 || round(xR1(i,5))==4
        for k=[-2 -1 0 1 2]
    II1{j+3,k+3,i}=ydatacrpdR1{i}(round(xR1(i,4))+k,round(xR1(i,2))+j);
        end
        end
    end
    end
    
        SII1=[II1{:,:,i}];
        xR1(i,6)=sum(SII1(:));
        
        
%---------- Third Spot ---------------------------------------------------
    if round(xR2(i,3))==2 || round(xR2(i,3))==1
        
        % shift x when it is too close to the left border
        if round(xR2(i,2))==0 || round(xR2(i,2))==1  
        xR2(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if round(xR2(i,2))==XSize(i) || round(xR2(i,2))==XSize(i)-1 
        xR2(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(xR2(i,5))==2 || round(xR2(i,5))==1
        for k=[-1 0 1]
    II2{j+2,k+2,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
        if round(xR2(i,5))==3 || round(xR2(i,5))==4
        for k=[-2 -1 0 1 2]
    II2{j+2,k+3,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
    end
    end
    
    if round(xR2(i,3))==3 || round(xR2(i,3))==4
        
        % shift x when it is too close to the left border
        if round(xR2(i,2))==0 || round(xR2(i,2))==1 || round(xR2(i,2))==2 
        xR2(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(xR2(i,2))==XSize(i) || round(xR2(i,2))==XSize(i)-1 || round(xR2(i,2))==XSize(i)-2 
        xR2(i,2)=XSize(i)-3;
        end
        
        
    for j=[-2 -1 0 1 2] 
        if round(xR2(i,5))==2 || round(xR2(i,5))==1
        for k=[-1 0 1]
    II2{j+3,k+2,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
        if round(xR2(i,5))==3 || round(xR2(i,5))==4
        for k=[-2 -1 0 1 2]
    II2{j+3,k+3,i}=ydatacrpdR2{i}(round(xR2(i,4))+k,round(xR2(i,2))+j);
        end
        end
    end
    end
    
        SII2=[II2{:,:,i}];
        xR2(i,6)=sum(SII2(:));
        
%---------- Fourth Spot ---------------------------------------------------

    if round(xR3(i,3))==2 || round(xR3(i,3))==1

        % shift x when it is too close to the left border
        if round(xR3(i,2))==0 || round(xR3(i,2))==1  
        xR3(i,2)=2;
        end
        
        % shift x when it is too close to the right border
        if round(xR3(i,2))==XSize(i) || round(xR3(i,2))==XSize(i)-1 
        xR3(i,2)=XSize(i)-2;
        end
        
        
    for j=[-1 0 1] 
        if round(xR3(i,5))==2 || round(xR3(i,5))==1
        for k=[-1 0 1]
    II3{j+2,k+2,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
        if round(xR3(i,5))==3 || round(xR3(i,5))==4
        for k=[-2 -1 0 1 2]
    II3{j+2,k+3,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
    end
    end
    
    if round(xR3(i,3))==3 || round(xR3(i,3))==4
    
        % shift x when it is too close to the left border
        if round(xR3(i,2))==0 || round(xR3(i,2))==1 || round(xR3(i,2))==2 
        xR3(i,2)=3;
        end
        
        % shift x when it is too close to the right border
        if round(xR3(i,2))==XSize(i) || round(xR3(i,2))==XSize(i)-1 || round(xR3(i,2))==XSize(i)-2 
        xR3(i,2)=XSize(i)-3;
        end
    
    for j=[-2 -1 0 1 2] 
        if round(xR3(i,5))==2 || round(xR3(i,5))==1
        for k=[-1 0 1]
    II3{j+3,k+2,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
        if round(xR3(i,5))==3 || round(xR3(i,5))==4
        for k=[-2 -1 0 1 2]
    II3{j+3,k+3,i}=ydatacrpdR3{i}(round(xR3(i,4))+k,round(xR3(i,2))+j);
        end
        end
    end
    end
    
        SII3=[II3{:,:,i}];
        xR3(i,6)=sum(SII3(:));
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