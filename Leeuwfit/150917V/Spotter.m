function [x0,Case,ydatacrpdR1,Ydata,Size,Yg] = Spotter(ydatacrpd,SA,Sx,Sy,Px,Py,Bs,lob,upb)
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
Case=2; 
elseif (Xg-SA)<=0
Ydata=ydatacrpd((Yg-SA):(Yg+SA),1:2*SA+1);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
Case=-1;
ydatacrpdR1((Yg-Bs):(Yg+Bs),1:Bss)=0;
elseif (Xg+SA)>XSize
Ydata=ydatacrpd((Yg-SA):(Yg+SA),XSize-2*SA:XSize);
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
ydatacrpdR1((Yg-Bs):(Yg+Bs),XSize-Bss+1:XSize)=0;
Case=1;
else 
Ydata=ydatacrpd((Yg-SA):(Yg+SA),(Xg-SA):(Xg+SA));
Ydata=[zeros(2*SA+1,Px) Ydata zeros(2*SA+1,Px)];
ydatacrpdR1((Yg-Bs):(Yg+Bs),(Xg-Bs):(Xg+Bs))=0;
Case=0;
end

[Dummy,I3]=max(Ydata);
[~,I4]=max(Dummy);
Xg2=I4;
Yg2=I3(I4);
x0=[Ampguess,Xg2,Sx,Yg2,Sy];% initial value


