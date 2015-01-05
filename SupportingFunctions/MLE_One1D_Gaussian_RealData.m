function paramsF=MLE_One1D_Gaussian_RealData(x,Y,est, plotit);
%

lx=length(x);

%%%This function fits one Gaussians to a line

% correctCoeffs = [ux1 ux2 sigmaValueX1 sigmaValueX2 b N1 N2];
a=1; % this should in principle be the pixel size
%Defining the initial conditions for the optimization procedure
ux50=x(ceil(0.5*lx));   %value x-axis on 0.2 of total length

Sx50=est.psf;  %pointspread in pixels

b0=0;          %background 
N50=max(Y)/(a^2*(1./(sqrt(2.*pi).*Sx10)));  %estimate based on half the signal max (assuming overlapping peaks)

x0 = [ux50,Sx50,b0,N50];

 % The function to be minimized is the negative of the log likelihood
datafun = @(params)((sum (expected (x,params,a))) - sum (Y.*log (expected (x,params,a))));
options = optimset ('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun', 1e-5);
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options); %we make use of fminsearch to find the minimum of the -log-likelihood function


%Here we generate data with the fitted parameters.
yFit = paramsF(4)*a^2*(1./(sqrt(2.*pi).*paramsF(2))).*exp((-((x-paramsF(1)).^2./(2.*paramsF(2).^2)))); %Generate y-values with the fitted parameters

 
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
    %close(gcf);
end
end


function p = oneDGauss (x,ux,s)
%This is the equation for a 1D gaussian
p = 1/(sqrt(2*pi).*s) .* exp (-((x-ux).^2./(2*s.^2) ));

end


function E = expected (x,params,a)
% This is the expected value given the parameters
E = abs(params (4).*a^2.*oneDGauss (x,params (1),params (2))) + params (3)^2 ;
end


