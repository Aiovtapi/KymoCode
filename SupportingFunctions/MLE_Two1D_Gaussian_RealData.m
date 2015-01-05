function paramsF=MLE_Two1D_Gaussian_RealData(x,Y,est, plotit);
%

lx=length(x);

%%%This function fits two Gaussians to a line

% correctCoeffs = [ux1 ux2 sigmaValueX1 sigmaValueX2 b N1 N2];
%Defining the initial conditions for the optimization procedure
ux10=est.x0;   %
ux20=est.x1;   %
b0=0;          %background 
Sx10=est.psf; Sx20=est.psf;
N10=est.N;  
N20=est.N;  
x0 = [ux10,ux20,Sx10,Sx20,b0,N10,N20];

 % The function to be minimized is the negative of the log likelihood
datafun = @(params)((sum (expected (x,params))) - sum (Y.*log (expected (x,params))));
options = optimset ('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun', 1e-5);
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options); %we make use of fminsearch to find the minimum of the -log-likelihood function


%Here we generate data with the fitted parameters.
y1Fit = paramsF(6).*exp((-((x-paramsF(1)).^2./(2.*paramsF(3).^2)))); %Generate y-values with the fitted parameters
y2Fit = paramsF(7).*exp((-((x-paramsF(2)).^2./(2.*paramsF(4).^2)))); %Generate y-values with the fitted parameters
YFit = y1Fit+y2Fit+b0; %add the generated y-values to obtain the final curve

 
if plotit==1
    %figure;
    plot(x,Y,'o-r','LineWidth',1.5);
    hold on;
    plot(x,y1Fit,'b-*','LineWidth',1.5);
    plot(x,y2Fit,'g-*','LineWidth',1.5);
    plot(x,(y1Fit+y2Fit),'k-*','LineWidth',1.5);
    % 
    h = legend('Input Data','Fit of Data 1','Fit of Data 2','Fit1+Fit2',5);
    set(h,'FontSize',14);
    xlabel ('Position ','fontsize',21);
    ylabel ('Intensity [a.u.]','fontsize',21);
    hold off;
    pause(1)
    [~]=ginput(1);
    %close(gcf);
end
end


function p = oneDGauss (x,ux,s)
%This is the equation for a 1D gaussian
p = exp (-((x-ux).^2./(2*s.^2) ));

end


function E = expected (x,params)
% This is the expected value given the parameters
E = (abs(params (6).*oneDGauss (x,params (1),params (3)))  +  abs(params (7).*oneDGauss (x,params(2),params(4))) + params (5)) ;
end


%from simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %clc;
% clear all;
% %format long;
% 
% 
% %Defining the 
% ux1=3.75; %mean of Gaussian 1
% ux2=5; %mean of Gaussian 2
% sigmaValueX1=0.5; %sigma of Gaussian 1
% sigmaValueX2=0.5; %sigma of Gaussian 1
% b=0; %this is in princple the background value
% N1=10000;  %the amplitude of Gaussian 1
% N2= 7000;; %the amplitude of Gaussian 2

% correctCoeffs = [ux1 ux2 sigmaValueX1 sigmaValueX2 b N1 N2];
% 
% %The x-values which we use to simulate the Gaussian with
% x=[1:0.1:9];
% 
% %Here we generate the data which we try to fit
% y1 = N1*a^2*(1./(sqrt(2.*pi).*sigmaValueX1)).*exp((-((x-ux1).^2./(2.*sigmaValueX1.^2)))); %Substituting the xValue in the 1D Gaussian function
% y2 = N2*a^2*(1./(sqrt(2.*pi).*sigmaValueX2)).*exp((-((x-ux2).^2./(2.*sigmaValueX2.^2)))); %Substituting the xValue in the 1D Gaussian function
% Y = y1+y2; % adding the two gaussians

