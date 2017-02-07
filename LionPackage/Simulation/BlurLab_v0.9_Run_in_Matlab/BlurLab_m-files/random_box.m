%Tristan Ursell
%Random 3D Diffusion of particles

function random_box(nframes,pts,meanI,D,Lx,Ly,Lz)
%create file name
tmp1=clock;
fname=['BlurLab_rand_output_' date '_' num2str(tmp1(4)) '-' num2str(tmp1(5)) '-' num2str(tmp1(6))...
    '_pts-' num2str(pts) '_meanI-' num2str(meanI) '_D-' num2str(D) '_Lx-' num2str(Lx) '_Ly-' num2str(Ly) '_Lz-' num2str(Lz) '.txt'];

nlines=nframes*pts;

X0=Lx*rand(1,pts);
Y0=Ly*rand(1,pts);
Z0=Lz*rand(1,pts);
L=1:pts;
I=poissrnd(meanI,1,pts);

Xout=zeros(1,nlines);
Yout=zeros(1,nlines);
Zout=zeros(1,nlines);
Iout=zeros(1,nlines);
Fout=zeros(1,nlines);
Lout=zeros(1,nlines);

for i=1:nframes    
    Xout((i-1)*pts+1:i*pts)=X0;
    Yout((i-1)*pts+1:i*pts)=Y0;
    Zout((i-1)*pts+1:i*pts)=Z0;
    Iout((i-1)*pts+1:i*pts)=I;
    Fout((i-1)*pts+1:i*pts)=i;
    Lout((i-1)*pts+1:i*pts)=L;
    
    phi=2*pi*rand(1,pts);
    theta=pi-acos(2*rand(1,pts)-1);
    R=normrnd(0,sqrt(2*D),1,pts);
    X0=X0+R.*sin(theta).*sin(phi);
    Y0=Y0+R.*sin(theta).*cos(phi);
    Z0=Z0+R.*cos(theta);
    
    %I know this isn't quite right
    X0(X0<0)=0;
    X0(X0>Lx)=Lx;
    Y0(Y0<0)=0;
    Y0(Y0>Ly)=Ly;
    Z0(Z0<0)=0;
    Z0(Z0>Lz)=Lz;
end

%write output file
BlurLab_text(Xout',Yout',Zout',Iout',Fout',Lout',fname)



