function MLE_Double2D_Gaussian_Final
clc;
clear all;
format long;

pxSize = 1; %um;
b=0;

ux1=pxSize*3;
uy1=pxSize*3;
sigmaValueX1=pxSize*1;
sigmaValueY1=pxSize*1;
N1=10^4;;
a1=1*pxSize;


xShift = 2;
yShift = 2;
ux2=pxSize*(3-xShift);
uy2=pxSize*(3-yShift);
sigmaValueX2=pxSize*1;
sigmaValueY2=pxSize*1;
N2=10^4;;
a1=1*pxSize;


correctCoeffs1 = [ux1 uy1 sigmaValueX1 sigmaValueY1,b,N1];
correctCoeffs2 = [ux2 uy2 sigmaValueX2 sigmaValueY2,b,N2];

%[x,y] = meshgrid(-9:0.25:9,-9:0.25:9);
[x,y] = meshgrid(-9*pxSize:pxSize/2:9*pxSize,-9*pxSize:pxSize/2:9*pxSize);
zValues1 = N1*a1^2*(1./(2.*pi.*sigmaValueX1.*sigmaValueY1)).*exp(-((x-ux1).^2./(2.*sigmaValueX1.^2) + (y-uy1).^2./(2.*sigmaValueY1.^2)));
zValues2 = N2*a1^2*(1./(2.*pi.*sigmaValueX2.*sigmaValueY2)).*exp(-((x-ux2).^2./(2.*sigmaValueX2.^2) + (y-uy2).^2./(2.*sigmaValueY2.^2)));

%corrupt with Poisson noise 
% J = imnoise(zValues,'gaussian',1,0.001)
% zValues=J;
% snr = 0.1;
% s_noisy1 = awgn(zValues1,snr);
% fn1 = zValues1 + randn(size(zValues1))*10; % Adding white noise
% 
% s_noisy2 = awgn(zValues2,snr);
% fn2 = zValues2 + randn(size(zValues2))*10; % Adding white noise

%fn = fn1 + fn2;

fn = zValues1 + zValues2;

figure
%mesh(x,y,zValues);
mesh(x,y,fn);


ux10=2;
uy10=2;
Sx10=0.5;
Sy10=0.5;
b0=0;
N10=10^4;

ux20=2;
uy20=2;
Sx20=0.5;
Sy20=0.5;
b0=0;
N20=10^4;


%x0 = [ux10, uy10, Sx10,Sy10,b0,N10];
x0 = [ux10, ux20, uy10, uy20, Sx10, Sx20, Sy10, Sy20, b0, N10, N20];

% The funtion to be minimized is the negative of the log likelihood
datafun = @(params)(sum (sum ((expected (x,y,params,a1))))-sum (sum (fn.*log (expected (x,y,params,a1)))));
options = optimset ('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun', 1e-5);
% fminsearch performs the multivariable minimization
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options);

correctCoeffs1
correctCoeffs2

paramsF


% zValues1Fit = paramsF (10)*a1^2*(1./(2.*pi.*paramsF (5).*paramsF (7))).*exp(-((x-paramsF (1)).^2./(2.*paramsF (5).^2) + (y-paramsF (3)).^2./(2.*paramsF (7).^2)));
% zValues2Fit = paramsF (11)*a1^2*(1./(2.*pi.*paramsF (6).*paramsF (8))).*exp(-((x-paramsF (2)).^2./(2.*paramsF (6).^2) + (y-paramsF (4)).^2./(2.*paramsF (8).^2)));

%YFit = zValues1Fit + zValues2Fit;                  
%YFit = paramsF (10)*a1^2*twoDGauss (x,y,paramsF (1),paramsF (3),paramsF (5),paramsF (7)) + paramsF (11)*a1^2*twoDGauss (x,y,paramsF (2),paramsF (4),paramsF (6),paramsF (8))+paramsF (9)^2;
% figure;
% mesh(x,y,YFit);
figure;
mesh (expected (x,y,paramsF,a1))
% paramsF = [3 3 3 3 1 1 1 1 0 10^4 10^4]
% mesh (expected (x,y,paramsF,a1))

end

function p = twoDGauss (x,y,ux,uy,s1,s2)
% Double 2D Gaussian
p = 1/(2.*pi.*s1.*s2) .* exp (-( (x-ux).^2./(2.*s1.^2)+(y-uy).^2./(2.*s2.^2) ));
end

function E = expected (x,y,params,a)
% The expected counts per pixel.
E = params (10).*a.^2.*twoDGauss (x,y,params (1),params (3),params (5),params (7)) + params (11).*a.^2.*twoDGauss (x,y,params (2),params (4),params (6),params (8)) + params (9).^2;
%E = params (10).*a.^2.*twoDGauss (x,y,params (1),params (3),params (5),params (7)) + params (11).*a.^2.*twoDGauss (x,y,params (2),params (4),params (6),params (8));
end


