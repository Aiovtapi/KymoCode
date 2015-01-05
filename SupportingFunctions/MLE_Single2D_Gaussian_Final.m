function MLE_Single2D_Gaussian_Final
clc;
clear all;
format long;

pxSize = 1; %um;

ux=pxSize*3;
uy=pxSize*3;
sigmaValueX=pxSize*1;
sigmaValueY=pxSize*1;
b=0;
N=10^4;;
a=1*pxSize;

correctCoeffs = [ux uy sigmaValueX sigmaValueY,b,N];

%[x,y] = meshgrid(-9:0.25:9,-9:0.25:9);
[x,y] = meshgrid(-9*pxSize:pxSize/2:9*pxSize,-9*pxSize:pxSize/2:9*pxSize);
zValues = N*a^2*(1./(2.*pi.*sigmaValueX.*sigmaValueY)).*exp(-((x-ux).^2./(2.*sigmaValueX.^2) + (y-uy).^2./(2.*sigmaValueY.^2)));


% %corrupt with Poisson noise 
% % J = imnoise(zValues,'gaussian',1,0.001)
% % zValues=J;
% snr = 0.1;
% s_noisy = awgn(zValues,snr);
% fn = zValues + randn(size(zValues))*10; % Adding white noise

fn = zValues;

figure
%mesh(x,y,zValues);
mesh(x,y,fn);

Sx= sigmaValueX;
Sy= sigmaValueY;

ux0=2;
uy0=2;
Sx0=0.5;
Sy0=0.5;
b0=0;
N0=50;


x0 = [ux0, uy0, Sx0,Sy0,b0,N0];

% The funtion to be minimized is the negative of the log likelihood
datafun = @(params)(sum (sum ((expected (x,y,params,a))))-sum (sum (zValues.*log (expected (x,y,params,a)))));
options = optimset ('MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-5);
% fminsearch performs the multivariable minimization
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options);

correctCoeffs

paramsF


figure;
mesh (expected (x,y,paramsF,a))

end

function p = twoDGauss (x,y,ux,uy,s1,s2)
% 2D Gaussian. (SOM Eq. 3)
p = 1/(2*pi*s1*s2) * exp (-( (x-ux).^2./(2*s1^2)+(y-uy).^2./(2*s2^2) ));

end

function E = expected (x,y,params,a)
% The expected counts per pixel. (SOM Eq. 12)
E = params (6)*a^2*twoDGauss (x,y,params (1),params (2),params (3),params (4))+params (5)^2;
end


