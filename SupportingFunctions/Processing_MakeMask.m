function mask=Processing_MakeMask(im,dif,psf);

[r,c]=size(im);
[x,y]=meshgrid(1:c,1:r);
xm=c/2;ym=r/2;
r0=((x-xm).^2+(y-ym).^2).^0.5; %'radial picture'
cutoff=c/4;
sa=abs(dif)+3*psf; sb=3*psf;
mask=exp(-(((x-xm)/sa).^4+((y-ym)/sb).^4));
