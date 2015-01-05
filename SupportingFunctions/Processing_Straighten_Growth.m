function [static_kymo,xmap]=Processing_Straighten_Growth(kymo,initval);
%This function removes the effect of (expected) bacterial growth
%from a plot, to ease visibility of division
%Create a 'growth line map'

[r,c]=size(kymo);
xmap=zeros(r,c);
[X,Y]=meshgrid(1:c,1:r);
x0=1:c;  % this is the initial position axis of bacterial material

for f=1:r
xmap(f,:)=x0*2^((f-1)/initval.estimateddoublingtime); %in frames); %These indicate the expected positions
end
static_kymo=interp2(kymo,xmap,Y);
% pcolor(static_kymo); shading flat; colormap hot
% [~]=ginput(1);
% close(gcf);

  