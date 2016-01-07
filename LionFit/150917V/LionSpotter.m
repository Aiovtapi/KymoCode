function [x0,Case,ydatacrpdR1,Ydata,Size,Yg,Xg] = LionSpotter(ydatacrpd,SA,Sx,Sy,Px,Py,Bs,lob,upb)
% Spotter: main spot detection + spot cut-out function

Bss=2*Bs+1;

[YSize,XSize]=size(ydatacrpd);
Size=[YSize,XSize];
[Amp,I]=max(ydatacrpd);
[Amp2,I2]=max(Amp);
Ampguess=max(Amp2);
Xg=I2;
Yg=I(I2);

ydatacrpdR1=ydatacrpd;

if Yg<lob || Yg>upb 
Ydata=zeros(2*SA+1,2*Px+2*SA+1);
ydatacrpdR1=0;
Case=10; 
elseif (Xg-SA)<=0 && (Yg-SA)>0 && (Yg+SA)<=YSize
Ydata=ydatacrpd((Yg-SA):(Yg+SA),1:2*SA+1);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
Case=2;
ydatacrpdR1((Yg-Bs):(Yg+Bs),1:Bss)=0;
elseif (Xg-SA)<=0 && (Yg+SA)>YSize
Ydata=ydatacrpd(YSize-2*SA:YSize,1:2*SA+1);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
Case=3;
ydatacrpdR1(YSize-Bss+1:YSize,1:Bss)=0;
elseif (Xg-SA)<=0 && (Yg-SA)<=0
Ydata=ydatacrpd(1:2*SA+1,1:2*SA+1);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];Case=-2;
ydatacrpdR1(1:Bss,1:Bss)=0;
Case=1;
elseif (Xg+SA)>XSize && (Yg-SA)>0 && (Yg+SA)<=YSize
Ydata=ydatacrpd((Yg-SA):(Yg+SA),XSize-2*SA:XSize);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
ydatacrpdR1((Yg-Bs):(Yg+Bs),XSize-Bss+1:XSize)=0;
Case=6;
elseif (Xg+SA)>XSize && (Yg+SA)>YSize
Ydata=ydatacrpd(YSize-2*SA:YSize,XSize-2*SA:XSize);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
ydatacrpdR1(YSize-Bss+1:YSize,YSize-Bss+1:XSize)=0;
Case=5;
elseif (Xg+SA)>XSize && (Yg-SA)<=0
Ydata=ydatacrpd(1:2*SA+1,XSize-2*SA:XSize);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
ydatacrpdR1(1:Bss,XSize-Bss+1:XSize)=0;
Case=7;
elseif (Yg+SA)>YSize && (Xg-SA)>0 && (Xg+SA)<=XSize
Ydata=ydatacrpd(YSize-2*SA:YSize,(Xg-SA):(Xg+SA));
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
ydatacrpdR1(YSize-Bss+1:YSize,(Xg-Bs):(Xg+Bs))=0;
Case=4;
elseif (Yg-SA)<=0 && (Xg-SA)>0 && (Xg+SA)<=XSize
Ydata=ydatacrpd(1:2*SA+1,(Xg-SA):(Xg+SA));
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
ydatacrpdR1(1:Bss,(Xg-Bs):(Xg+Bs))=0;
Case=8;
else 
Ydata=ydatacrpd((Yg-SA):(Yg+SA),(Xg-SA):(Xg+SA));
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Ydata=[zeros(Py,size(Ydata,2));Ydata;zeros(Py,size(Ydata,2))];
ydatacrpdR1((Yg-Bs):(Yg+Bs),(Xg-Bs):(Xg+Bs))=0;
Case=9;
end

[Dummy,I3]=max(Ydata);
[~,I4]=max(Dummy);
Xg2=I4;
Yg2=I3(I4);
x0=[Ampguess,Xg2,Sx,Yg2,Sy];% initial value


