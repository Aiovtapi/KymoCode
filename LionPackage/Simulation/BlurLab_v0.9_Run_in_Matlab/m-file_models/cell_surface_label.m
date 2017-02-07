%Tristan Ursell
%Cell Surface Diffusion

%last modified August 9th 2011

close all
clearvars
pause(0.01)

%Number of frames
M=100;

%Cell Radius (um)
R=0.5;

%Cell length (um)
L=2;

%Areal density of markers (1/um^2)
%max is 10000 to "fit"
rho=200;

%cap area (um^2)
Acap=2*pi*R^2;

%cylinder area (um^2)
Acyl=2*pi*R*L;

%Position sigma (um)
sig=0.02;

%total number
Nlhs=poissrnd(Acap*rho);
Nrhs=poissrnd(Acap*rho);
Ncyl=poissrnd(Acyl*rho);    

%initialize vectors
Xtot=[];
Ytot=[];
Ztot=[];
Itot=[];
Ftot=[];
Ltot=[];

for i=1:M
    
    if i==1
        %lhs hemispherical cap
        theta1=pi-acos(2*rand(1,Nlhs)-1);
        phi1=pi+pi*rand(1,Nlhs);
        X1=R*sin(theta1).*sin(phi1);
        Y1=R*sin(theta1).*cos(phi1);
        Z1=R*cos(theta1);
        
        %rhs hemispherical cap
        theta2=pi-acos(2*rand(1,Nrhs)-1);
        phi2=pi*rand(1,Nrhs);
        X2=L+R*sin(theta2).*sin(phi2);
        Y2=R*sin(theta2).*cos(phi2);
        Z2=R*cos(theta2);
        
        %cylindrical section
        theta3=2*pi*rand(1,Ncyl);
        X3=L*rand(1,Ncyl);
        Y3=R*cos(theta3);
        Z3=R*sin(theta3);
    end
    
    %random motion (constrained)
    Ntot=Nlhs+Nrhs+Ncyl;
    Rnd=normrnd(0,sig,1,Ntot);
    thetand=pi-acos(2*rand(1,Ntot)-1);
    phind=2*pi*rand(1,Ntot);
    Xnd=Rnd.*sin(thetand).*sin(phind);
    Ynd=Rnd.*sin(thetand).*cos(phind);
    Znd=Rnd.*cos(thetand);
    
    %Current positions
    Xcurr=[X1,X2,X3]+Xnd;
    Ycurr=[Y1,Y2,Y3]+Ynd;
    Zcurr=[Z1,Z2,Z3]+Znd;
    
    Ntemp=length([X1,X2,X3]);
    Xtot=[Xtot,Xcurr];
    Ytot=[Ytot,Ycurr];
    Ztot=[Ztot,Zcurr];
    Itot=[Itot,ones(1,Ntot)];
    Ftot=[Ftot,i*ones(1,Ntot)];
    Ltot=[Ltot,1:Ntot];
    
    disp(num2str(i))
end

%generate a cell
Npts=200;
Reff=R-0.01;
[Xsph,Ysph,Zsph]=sphere(Npts);
[Zcyl,Ycyl,Xcyl]=cylinder(Reff,2*Npts);

figure;
mksz=15;
hold on
surf(L*Xcyl,Ycyl,Zcyl)
surf(Reff*Xsph,Reff*Ysph,Reff*Zsph)
surf(Reff*Xsph+L,Reff*Ysph,Reff*Zsph)
shading('interp')
colormap(copper)
colormap([1,0.6324,0.4027])
%plot3(Xcurr,Ycurr,Zcurr,'.','MarkerSize',mksz)
plot3([X1,X2,X3],[Y1,Y2,Y3],[Z1,Z2,Z3],'.','MarkerSize',mksz,'color',[0.1 0.8 0.1])
axis image
axis vis3d
box off
view([-50,25])
lightangle(50,60)
lighting gouraud
material default
xlim([-R-0.1,L+R+0.1])
ylim([-R-0.1,R+0.1])
zlim([-R-0.1,R+0.1])
xlabel('X')
ylabel('Y')
zlabel('Z')

%creat mesh
cvec=[0 0 0];
nlines=10;
zline=0;
ylines=linspace(-R-0.1,R+0.1,nlines);
xlines=linspace(-R-0.1,L+R+0.1,round(nlines*(2*R+0.2+L)/(2*R+0.2)));
for i=1:length(xlines)
    plot3([xlines(i),xlines(i)],[ylines(1),ylines(end)],[zline,zline],'-','color',cvec,'LineWidth',2)
end
for i=1:length(ylines)
    plot3([xlines(1),xlines(end)],[ylines(i),ylines(i)],[zline,zline],'-','color',cvec,'LineWidth',2)
end

%creat mesh
nlines=10;
zline=0.4;
ylines=linspace(-R-0.1,R+0.1,nlines);
xlines=linspace(-R-0.1,L+R+0.1,round(nlines*(2*R+0.2+L)/(2*R+0.2)));
for i=1:length(xlines)
    plot3([xlines(i),xlines(i)],[ylines(1),ylines(end)],[zline,zline],'-','color',cvec,'LineWidth',2)
end
for i=1:length(ylines)
    plot3([xlines(1),xlines(end)],[ylines(i),ylines(i)],[zline,zline],'-','color',cvec,'LineWidth',2)
end

%figure;
%generate a cell
Npts=20;
[Xsph,Ysph,Zsph]=sphere(Npts);
[Zcyl,Ycyl,Xcyl]=cylinder(R,2*Npts);

%surf(0.1*Xcyl,Ycyl,Zcyl)
%surf(R*Xsph,R*Ysph,R*Zsph)
%axis image
%axis vis3d
%box on
%view([-50,25])

q1=input('Create file? ','s');

if strcmp(q1,'yes')
    fname=[date '_cell_surface_R-' num2str(R) '_rho-' num2str(rho) '.txt'];
    imsim_text(Xtot',Ytot',Ztot',Itot',Ftot',Ltot',fname)
    disp(['Wrote file with name: ' fname])
end

%*****************************************************
break

[file1,aa]=imgetfile;

Nfile=length(imfinfo(file1));

Im=zeros([size(imread(file1,1)),Nfile]);

for i=1:Nfile
    Im(:,:,i)=imread(file1,i);
end

Imean=mean(Im,3);
Istd=std(Im,0,3);

%for 2000
%cutoff=2200;
%frac=0.01;

%for 500
%cutoff=550;
%frac=0.05;

%for 100
%cutoff=110;
%frac=0.1;

%for top 500
cutoff=500;
frac=0.1;

Ibase=Imean;
Ibase(Ibase<=cutoff)=cutoff;

cig=60;
rm=200;
gm=100;
bm=30;
tsumap=zeros(256,3);
vals=0:1:255;
tsumap(:,1)=exp(-(vals-rm).^2/(2*cig^2));
tsumap(:,2)=exp(-(vals-gm).^2/(3*cig^2));
tsumap(:,3)=exp(-(vals-bm).^2/(3*cig^2));

h1=figure;
axes1=axes('Parent',h1,'ZMinorTick','on','DataAspectRatio',[1 1 1/frac]);
hold on
shading interp
surf([Ibase(15:35,10:60);Istd(15:35,10:60)+cutoff]-cutoff,'Parent',axes1)
axis tight
box on
light
lightangle(0,60)
lighting gouraud
view([-120,25])
view([-80,25])
xlabel('X')
ylabel('Y')
zlabel('Int')
colormap(tsumap)

%****************************************************
%make nice DUAL airy disk
lambda=0.532;  %nm
k=2*pi/lambda;

%aperature
a=0.1;

%distance from source
d=10;

%separation (um)
sep=80*lambda*0;

%step size
dx=round(1/100*d^2/lambda);

%square limits
lim=12*d*lambda;

%initialize vectors
Xvec=(-lim-sep/2):dx:(lim+sep/2);
Yvec=(-lim):dx:(lim);

%setup matrics
Xairy=zeros(length(Yvec),length(Xvec));
Yairy=zeros(length(Yvec),length(Xvec));
for i=1:length(Xvec)
    Yairy(:,i)=Yvec;
end
for i=1:length(Yvec)
    Xairy(i,:)=Xvec;
end

%calculate Airy
Zairy=(besselj(1,k*a/d*sqrt(1e-12+(Xairy-sep/2).^2+Yairy.^2))./(k*a/d*sqrt(1e-12+(Xairy-sep/2).^2+Yairy.^2))).^2+...
    (besselj(1,k*a/d*sqrt(1e-12+(Xairy+sep/2).^2+Yairy.^2))./(k*a/d*sqrt(1e-12+(Xairy+sep/2).^2+Yairy.^2))).^2;

sig=0.01;
Znoise=0*normrnd(4*sig,sig,size(Zairy));

h2=figure;
axes1=axes('Parent',h2,'ZMinorTick','on','DataAspectRatio',[1 1 0.005]);
hold on
shading interp
surf(Xairy,Yairy,Zairy+Znoise,'Parent',axes1)
axis tight
box on
light
lightangle(0,0)
lighting gouraud
view([-30,20])
xlabel('X')
ylabel('Y')
zlabel('Int')
colormap(tsumap)


Zgauss1=exp(-(Xairy.^2+Yairy.^2)/(2*4^2));
Zgauss2=exp(-(Xairy.^2+Yairy.^2)/(2*20^2));
[Zg1x,Zg1y]=gradient(Zgauss1);
[Zg2x,Zg2y]=gradient(Zgauss2);

Z1mag=sqrt(Zg1x.^2+Zg1y.^2);
Z2mag=sqrt(Zg2x.^2+Zg2y.^2);

zmin=0;
zmax=max(max(Z2mag));

figure;
colormap(gray)
subplot(2,2,1)
imagesc(Zgauss2)
axis equal
axis tight

subplot(2,2,2)
imagesc(Zgauss1)
axis equal
axis tight

subplot(2,2,3)
imagesc(Z2mag,[zmin,zmax])
axis equal
axis tight

subplot(2,2,4)
imagesc(Z1mag,[zmin,zmax])
axis equal
axis tight










