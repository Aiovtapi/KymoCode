function MLE_Two1D_Gaussian_1Ver
clc;
clear all;
%format long;


%Defining the 
ux1=3.75; %mean of Gaussian 1
ux2=5; %mean of Gaussian 2
sigmaValueX1=0.5; %sigma of Gaussian 1
sigmaValueX2=0.5; %sigma of Gaussian 1
b=0; %this is in princple the background value
N1=10000;  %the amplitude of Gaussian 1
N2= 7000;; %the amplitude of Gaussian 2
a=1; % this should in principle be the pixel size
correctCoeffs = [ux1 ux2 sigmaValueX1 sigmaValueX2 b N1 N2];

%The x-values which we use to simulate the Gaussian with
x=[1:0.1:9];

%Here we generate the data which we try to fit
y1 = N1*a^2*(1./(sqrt(2.*pi).*sigmaValueX1)).*exp((-((x-ux1).^2./(2.*sigmaValueX1.^2)))); %Substituting the xValue in the 1D Gaussian function
y2 = N2*a^2*(1./(sqrt(2.*pi).*sigmaValueX2)).*exp((-((x-ux2).^2./(2.*sigmaValueX2.^2)))); %Substituting the xValue in the 1D Gaussian function
Y = y1+y2; % adding the two gaussians


%Defining the initial conditions for the optimization procedure
ux10=2;
ux20=8;
Sx10=0.5;
Sx20=0.5;
b0=0;
N10=9000;
N20=6000;
x0 = [ux10,ux20,Sx10,Sx20,b,N10,N20];

 % The funtion to be minimized is the negative of the log likelihood
datafun = @(params)((sum (expected (x,params,a))) - sum (Y.*log (expected (x,params,a))));
options = optimset ('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun', 1e-5);
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options); %we make use of fminsearch to find the minimum of the -log-likelihood function

correctCoeffs

paramsF

%Here we generate data with the fitted parameters.
y1Fit = paramsF(6)*a^2*(1./(sqrt(2.*pi).*paramsF(3))).*exp((-((x-paramsF(1)).^2./(2.*paramsF(3).^2)))); %Generate y-values with the fitted parameters
y2Fit = paramsF(7)*a^2*(1./(sqrt(2.*pi).*paramsF(4))).*exp((-((x-paramsF(2)).^2./(2.*paramsF(4).^2)))); %Generate y-values with the fitted parameters
YFit = y1Fit+y2Fit; %add the generated y-values to obtain the final curve

 
figure;
plot(x,Y,'o-r','LineWidth',1.5);
hold on;
plot(x,y1Fit,'b*','LineWidth',1.5);
plot(x,y2Fit,'g*','LineWidth',1.5);
plot(x,(y1Fit+y2Fit),'k*','LineWidth',1.5);
% 
h = legend('Simulated Data','Fit of Data 1','Fit of Data 2','Fit1+Fit2',5);
set(h,'FontSize',14);
xlabel ('Position ','fontsize',21);
ylabel ('Intensity [a.u.]','fontsize',21);
hold off;

end


function p = oneDGauss (x,ux,s)
%This is the equation for a 1D gaussian
p = 1/(sqrt(2*pi).*s) .* exp (-((x-ux).^2./(2*s.^2) ));

end


function E = expected (x,params,a)
% This is the expected value given the parameters
E = (params (6).*a^2.*oneDGauss (x,params (1),params (3)))  +  (params (7).*a^2.*oneDGauss (x,params(2),params(4))) + params (5)^2 ;
end
