clear all
close all

imSize = 300;                           % image size: n X n
sigma = 10;
gauss=cell(1,10);

% make linear ramp
X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
[Xm,Ym] = meshgrid(X0, X0);
s = sigma / imSize;  

for i=1:10
gauss{i}= i*exp( -((((Xm+(i-5)/10).^2)+((Ym).^2)) ./ (2* s^2)) ); % formula for 2D gaussian
end
Y=gauss{4};
X=imagesc(gauss{4});                        % display
colormap gray(256);
axis off; axis image;     % use gray colormap

imwrite(Y,'test','TIFF');