%Tristan Ursell
%Spherical Diffusion Simulator for BlurLab

%last modified August 4 2011, 9pm

close all
clearvars
pause(0.01)

%Number of iterations
M=500;

%Sphere radius (um)
R=2;

%Areal density of markers (1/um^2)
rho=200;

%Sphere area (um^2)
A=4*pi*R^2;

%Diffusion coefficient (um^2/s)
D=1;

%time step (s)
dt=0.02;

%mean jump distance (flat space) (um)
dx=sqrt(4*D*dt);

%total number of particles
N=poissrnd(A*rho);

%initialize vectors
Xtot=zeros(1,N*M);
Ytot=zeros(1,N*M);
Ztot=zeros(1,N*M);
Itot=ones(1,N*M);
Ftot=zeros(1,N*M);
Ltot=zeros(1,N*M);

%initial spherical positions
theta=pi-acos(2*rand(1,N)-1);
phi=2*pi*rand(1,N);

%simulate spherical diffusion
%this algorithm is approximately correct for dx/R<<1
for i=1:M
    %convert points from spherical
    X=R*sin(theta).*sin(phi);
    Y=R*sin(theta).*cos(phi);
    Z=R*cos(theta);
    
    %constrained, random surface jitter
    Rnd=normrnd(0,dx,1,N);
    thetand=pi-acos(2*rand(1,N)-1);
    phind=2*pi*rand(1,N);
    Xnd=X+Rnd.*sin(thetand).*sin(phind);
    Ynd=Y+Rnd.*sin(thetand).*cos(phind);
    Znd=Z+Rnd.*cos(thetand);
    
    Rtmp=sqrt(Xnd.^2+Ynd.^2+Znd.^2);
    theta=acos(Znd./Rtmp);
    phi=atan2(Xnd,Ynd);
    
    %create output vectors
    Xtot((i-1)*N+1:i*N)=X;
    Ytot((i-1)*N+1:i*N)=Y;
    Ztot((i-1)*N+1:i*N)=Z;
    Ftot((i-1)*N+1:i*N)=i*ones(1,N);
    Ltot((i-1)*N+1:i*N)=1:N;
    
    plot3(X,Y,Z,'r.')
    axis equal
    pause(0.1)
    
    disp(num2str(i))
end

%write output file
q1=input('Create file? ','s');

if strcmp(q1,'yes')
    fname=[date '_spherical_diffusion_R-' num2str(R) '_rho-' num2str(rho) '_dx-' num2str(dx) '.txt'];
    BlurLab_text(Xtot',Ytot',Ztot',Itot',Ftot',Ltot',fname)
    disp(['Wrote file with name: ' fname])
end



