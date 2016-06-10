%Tristan Ursell
%Cell Surface Marker for ImSim

%last modified June 1 2011, 6pm

close all
clearvars
pause(0.01)

%Number of iterations
M=100;

%Cell Radius (um)
R=0.5;

%Cell length (um)
L=2;

%Areal density of markers (1/um^2)
%max is 10000 to "fit"
rho=500;

%cap area (um^2)
Acap=2*pi*R^2;

%cylinder area (um^2)
Acyl=2*pi*R*L;

%Position sigma (um)
sig=0.02;

%initialize vectors
Xtot=[];
Ytot=[];
Ztot=[];
Itot=[];
Ftot=[];

for i=1:M
    
    %lhs hemispherical cap
    Nlhs=poissrnd(Acap*rho);
    theta1=pi-acos(2*rand(1,Nlhs)-1);
    phi1=pi+pi*rand(1,Nlhs);
    X1=R*sin(theta1).*sin(phi1);
    Y1=R*sin(theta1).*cos(phi1);
    Z1=R*cos(theta1);
    
    %rhs hemispherical cap
    Nrhs=poissrnd(Acap*rho);
    theta2=pi-acos(2*rand(1,Nrhs)-1);
    phi2=pi*rand(1,Nrhs);
    X2=L+R*sin(theta2).*sin(phi2);
    Y2=R*sin(theta2).*cos(phi2);
    Z2=R*cos(theta2);
    
    %cylindrical section
    Ncyl=poissrnd(Acyl*rho);
    theta3=2*pi*rand(1,Ncyl);
    X3=L*rand(1,Ncyl);
    Y3=R*cos(theta3);
    Z3=R*sin(theta3);
    
    %random motion
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
    
    disp(num2str(i))
end

%generate a cell
Npts=200;
Reff=R-0.01;
[Xsph,Ysph,Zsph]=sphere(Npts);
[Zcyl,Ycyl,Xcyl]=cylinder(Reff,2*Npts);

h1=figure;
axes1=axes('Parent',h1,'Fontsize',14,'ZTick',[-0.5 0 0.5],'YTick',[-0.5 0 0.5],'XTick',[0 1 2]);
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

% %creat mesh
% nlines=10;
% zline=0.4;
% ylines=linspace(-R-0.1,R+0.1,nlines);
% xlines=linspace(-R-0.1,L+R+0.1,round(nlines*(2*R+0.2+L)/(2*R+0.2)));
% for i=1:length(xlines)
%     plot3([xlines(i),xlines(i)],[ylines(1),ylines(end)],[zline,zline],'-','color',cvec,'LineWidth',2)
% end
% for i=1:length(ylines)
%     plot3([xlines(1),xlines(end)],[ylines(i),ylines(i)],[zline,zline],'-','color',cvec,'LineWidth',2)
% end


figure;
%generate a cell
Npts=20;
[Xsph,Ysph,Zsph]=sphere(Npts);
[Zcyl,Ycyl,Xcyl]=cylinder(R,2*Npts);

surf(0.1*Xcyl,Ycyl,Zcyl)
%surf(R*Xsph,R*Ysph,R*Zsph)
axis image
axis vis3d
box on
view([-50,25])

q1=input('Create file? ','s');

if strcmp(q1,'yes')
    fname=[date '_cell_surface_R-' num2str(R) '_rho-' num2str(rho) '.txt'];
    imsim_text(Xtot',Ytot',Ztot',Itot',Ftot',fname)
    disp(['Wrote file with name: ' fname])
end


