%Simulated Helix Maker
%Tristan Ursell
%May 2011

clear all
close all
clc
pause(0.1)

%Number of frames
N=100;

%Helix length (nm) (x-projected)
L=2000;

%Helix Y radius (nm)
Ry=500;

%Helix Z radius (nm)
Rz=Ry;

%Helix Pitch (turns/nm)
ptch=0.002;

%Stochasic density of emitters (1/nm arc)
lam=0.025;

%Helix rotation (rads)
phi=0;

%per frame diffusion coefficient for points (nm)
sig=10;

%rotation per frame (rads)
dtheta=2*pi/N;

%*********************************************************************

%Calculate rotation marix
Rot1=zeros(3,3);
Rot1(1,1)=cos(phi);
Rot1(1,3)=sin(phi);
Rot1(2,2)=1;
Rot1(3,1)=-sin(phi);
Rot1(3,3)=cos(phi);

%Emitter intensity (au)
intemt=1;

%Helix spatial resolution (nm)
dx=0.1;

%create output file name
fname=['helix_input_file_' date '.txt'];

%number of points in helix
Npts=length(0:dx:L);

if length(dir(fname))>0
    delete(fname)
end

figure;
fid=fopen(fname,'w+');
for j=1:N
    %make helix shape
    X0=0:dx:L;
    Y0=Ry*cos(2*pi*X0*ptch+(j-1)*dtheta);
    Z0=Rz*sin(2*pi*X0*ptch+(j-1)*dtheta);
    
    X=X0+normrnd(0,sig,1,Npts);
    Y=Y0+normrnd(0,sig,1,Npts);
    Z=Z0+normrnd(0,sig,1,Npts);
    
    %approximate helix point spacing
    ds=sqrt((X0(2)-X0(1))^2+(Y0(2)-Y0(1))^2+(Z0(2)-Z0(1))^2);
    
    if j==1
        %Pick emitter positions
        E=poissrnd(lam*ds,1,length(X))>0;
        
        %Create emitter vector
        int1=intemt*ones(1,sum(E));
    end
    
    for i=1:length(X)
        temp0=Rot1*[X0(i),Y0(i),Z0(i)]';
        Xrot0(i)=temp0(1);
        Yrot0(i)=temp0(2);
        Zrot0(i)=temp0(3);
        
        temp1=Rot1*[X(i),Y(i),Z(i)]';
        Xrot1(i)=temp1(1);
        Yrot1(i)=temp1(2);
        Zrot1(i)=temp1(3);
    end
    
    plot3(Xrot0,Yrot0,Zrot0,'k')
    hold on
    plot3(Xrot1(E),Yrot1(E),Zrot1(E),'.','color',[0 0.7 0],'Markersize',20)
    hold off
    axis equal
    axis tight
    pause(0.1)
    
    %write output file
    frm1=j*ones(1,sum(E));
    output=[Xrot1(E);Yrot1(E);Zrot1(E);int1;frm1];
    fprintf(fid,'%6.1f %6.1f %6.1f %6.1f %6f\n',output);
end
fclose(fid);



