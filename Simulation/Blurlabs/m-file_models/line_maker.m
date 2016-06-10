%Simulated Line Maker
%Tristan Ursell
%May 2011

clear all
close all
clc
pause(0.1)

%Number of frames
N=500;

%line length (nm) (x-projected)
L=3000;

%Stochasic density of emitters (1/nm arc)
lam=0.005;

%line pitch 
m=pi/2;

%line offset (nm)
b=100;

%line rotation (rads)
phi=pi/5;

%per frame diffusion coefficient for points (nm)
sig=0;

%speed of line translation (nm/frame)
spd=1;

%************************************************************************

%Calculate rotation marix
Rot1=zeros(3,3);
Rot1(1,1)=cos(phi);
Rot1(1,3)=sin(phi);
Rot1(2,2)=1;
Rot1(3,1)=-sin(phi);
Rot1(3,3)=cos(phi);

%Emitter intensity (au)
intemt=1;

%linear spatial resolution (nm)
dx=0.1;

%create output file name
fname=['line_input_file_' date '_L-' num2str(L) '_lam-' num2str(lam) '.txt'];

%number of points in helix
Npts=length(0:dx:L);

if ~isempty(dir(fname))
    fclose all;
    delete(fname)
end

figure;
fid=fopen(fname,'w+');
for j=1:N
    %make line
    X0=(0:dx:L)+(j-1)*spd;
    Y0=m*X0+b;
    Z0=zeros(1,Npts);
    
    X=X0+normrnd(0,sig,1,Npts);
    Y=Y0+normrnd(0,sig,1,Npts);
    Z=Z0+normrnd(0,sig,1,Npts);
    
    if j==1
        %approximate linear point spacing
        ds=sqrt((X0(2)-X0(1))^2+(Y0(2)-Y0(1))^2+(Z0(2)-Z0(1))^2);
        
        %Pick emitter positions
        E=poissrnd(lam*ds,1,length(X))>0;
        
        %Create emitter vector
        int1=intemt*ones(1,sum(E));
    end
    
    for i=1:length(X)
        temp0=Rot1*[X0(i),Y0(i),Z0(i)]';
        Xrot0(i)=temp0(1);
        Yrot0(i)=temp0(2);
        Zrot0(i)=temp0(3)+sin(phi)*L/2;
        
        temp1=Rot1*[X(i),Y(i),Z(i)]';
        Xrot1(i)=temp1(1);
        Yrot1(i)=temp1(2);
        Zrot1(i)=temp1(3)+sin(phi)*L/2;
    end
    
    plot3(Xrot0,Yrot0,Zrot0,'k')
    hold on
    plot3(Xrot1(E),Yrot1(E),Zrot1(E),'.','color',[0 0.7 0],'Markersize',20)
    hold off
    axis equal
    axis tight
    title([num2str(j) '/' num2str(N)])
    pause(0.001)
    
    %write output file
    %frm1=j*ones(1,sum(E));
    %output=[Xrot1(E);Yrot1(E);Zrot1(E);int1;frm1];
    %fprintf(fid,'%6.1f %6.1f %6.1f %6.1f %6f\n',output);
end
fclose(fid);



