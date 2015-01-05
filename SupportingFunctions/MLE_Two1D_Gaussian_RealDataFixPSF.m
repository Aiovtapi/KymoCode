function [paramsF,y1Fit,y2Fit,yAllFit]=MLE_Two1D_Gaussian_RealDataFixPSF(x,Y,est);
%%%This function fits two Gaussians to a line

%Defining the initial conditions for the optimization procedure
x0 = [est.x0,est.x1,est.N,est.N];

 % The function to be minimized is the negative of the log likelihood
datafun = @(params)((sum (expectedFixPSF (x,params,est.psf))) - sum (Y.*log (expectedFixPSF (x,params,est.psf))));
options = optimset ('MaxFunEvals', 1E3, 'MaxIter', 1E4, 'TolFun', 1e-4);
paramsF = fminsearch (datafun, x0, options); %we make use of fminsearch to find the minimum of the -log-likelihood function

paramsF(3)=abs(paramsF(3));  %to ensure positive amplitudes (also in expectation expression)
paramsF(4)=abs(paramsF(4));  %to ensure positive amplitudes (also in expectation expression)


%Here we generate data with the fitted parameters.
y1Fit =1/(est.psf*(2*pi)^0.5)*paramsF(3)*exp((-((x-paramsF(1)).^2./(2.*est.psf.^2)))); %Generate y-values with the fitted parameters
y2Fit = 1/(est.psf*(2*pi)^0.5)*paramsF(4)*exp((-((x-paramsF(2)).^2./(2.*est.psf.^2)))); %Generate y-values with the fitted parameters
yAllFit = y1Fit + y2Fit; %add the generated y-values to obtain the final curve

end


function p = oneDGauss (x,ux,s)
%This is the equation for a normalized1D gaussian
p = 1/(s*sqrt(2*pi))*exp (-(x-ux).^2./(2*s.^2));
end


function E = expectedFixPSF (x,params,psf)
% This is the expected value given the parameters
E = abs(params(3)).*oneDGauss (x,params (1),psf)  +  abs(params (4)).*oneDGauss (x,params(2),psf) ;
end

