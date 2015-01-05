function spotim=B201_Bacsimfocus(FL, focus, x, y, initval,plotit);
% This script builds a fake bacterium
% Jacob Kers 2012

if nargin <3  %TEST
 close all  
%We Start with ten-fold higher sampling (unit is 100 nm) (Ref Cheezum et
%al)
plotit=0;
FL=ones(20,32);
focus=1000;
initval.psf=1.3;
x=10;
y=10;
end

[r0,c0]=size(FL);
blowup=11;
r=r0*blowup;
c=c0*blowup;

hpsf=initval.psf*blowup;  %Point spread function
squ=1*blowup*hpsf;

xbl=x*blowup;  %two-third to one side
ybl=y*blowup;

[x,y]=meshgrid(1:c,1:r);

%2) Add a single-pixel emitter-------------------------------------
Bacpic=zeros*x;
x0=round(xbl);
y0=round(ybl);
Bacpic(y0,x0)=Bacpic(y0,x0)+focus;


%blur this picture with a typical kernel----------------------
[xk,yk]=meshgrid(1:squ,1:squ);
radius=((xk-squ/2).^2+(yk-squ/2).^2).^0.5; %'radial picture'
blurkernel=1/(2*pi*(hpsf.^2))*exp(-radius.^2/(2*hpsf.^2));
blurred_pic=imfilter(Bacpic,blurkernel);


%Throw this picture on a camera pixel grid
r0=round(r0);
c0=round(c0);
camera_pic=0*FL;
for i=1:r0;
    for j=1:c0
        loc=1+blowup*(j-1);
        hic=blowup*(j);
        lor=1+blowup*(i-1);
        hir=blowup*(i);
        squ=blurred_pic(lor:hir,loc:hic);
        camera_pic(i,j)=sum(squ(:));
    end
end


if plotit
% subplot(2,2,3); 
pcolor(camera_pic); shading flat, colormap hot;axis equal; title('downsampled');
sum(camera_pic(:))
[~]=ginput(1);
end
spotim=camera_pic;