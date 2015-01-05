function geom=Build_2D_Estimate(pic,x1,x2)
          


%-------------------------------Build 2D fit Estimate
  %find size pic
  
%   get crosssection per spot
%   Perform 1D gaussfit for y pos
%   now x1,y1
%   fix psf
%   get background
  
 %2x 1D MaxFit%%%%%%%%%%%%%%%%%%%%%5
 pic=pic-min(pic(:));
 pic=pic/max(pic(:));
 [r,c]=size(pic);
 x1r=ceil(x1); x1r=max([x1r 1]); x1r=min([x1r c]);  %check in-range
 x2r=ceil(x2); x2r=max([x2r 1]); x2r=min([x2r c]);
 yprof1=pic(:,x1r); [N1,y1]=max(yprof1);
 yprof2=pic(:,x2r); [N2,y2]=max(yprof2);
 

 b0=0;
 psf=1.2;
%  %%----------------------------------
%  %%analyze the image pattern
% [xm,ym,theta,ecc]=Processing_Calculate2Dmoment(pic.^4);
% 
% geom.theta=theta;
% geom.eccent=ecc;
%  
geom=1;
%  %-----------------------------
%  %I Do 2D fit using a double Gauss with open
%  %psfx and psfy 
%  % The function to be minimized is the negative of the log
%  % likelihood;  % fminsearch performs the multivariable minimization
% [x,y] = meshgrid(1:1:r,1:1:c);
% x0 = [y1, y2, x1, x2, psf, psf, psf, psf, b0, N1, N2];
% datafun = @(params)(sum (sum ((DoubleTwoDGaussexpected (x,y,params))))-sum (sum (pic'.*log (DoubleTwoDGaussexpected(x,y,params))))); 
% %I flip the image here "fn'" so as to have flipped coordinates at the top
% options = optimset ('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun', 1e-5);
% 
% [paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options);
% 
% frFit = DoubleTwoDGaussexpected(x,y,paramsF)';
% %--------------------------------------------------------

%-----------------------------
 %I Do 2D fit using a double Gauss with open
 %psfx and psfy 
 % The function to be minimized is the negative of the log
 % likelihood;  % fminsearch performs the multivariable minimization
[x,y] = meshgrid(1:1:r,1:1:c);
x0 = [y1, y2, x1, x2, b0, N1, N2];
datafun = @(params)(sum (sum ((DoubleTwoDGaussexpected_Fix1PSF(x,y,psf,params))))-sum (sum (pic'.*log (DoubleTwoDGaussexpected_Fix1PSF(x,y,psf,params))))); 
%I flip the image here "fn'" so as to have flipped coordinates at the top
options = optimset ('MaxFunEvals', 100000, 'MaxIter', 1000, 'TolFun', 1e-5);

[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options);

frFit = DoubleTwoDGaussexpected_Fix1PSF(x,y,psf,paramsF)';
x3=paramsF(3); y3=paramsF(1);
x4=paramsF(4); y4=paramsF(2);
%--------------------------------------------------------

sho=1;
if sho
%figure;
subplot(1,2,1);
 pcolor(pic); shading flat; colormap hot; title('2x1D Fit'); axis square ; hold on;
 plot(x3+0.5,y3+0.5,'o','MarkerFaceColor', 'b', 'MarkerSize', 8, 'MarkerEdgeColor','b');
 plot(x4+0.5,y4+0.5,'ro','MarkerFaceColor', 'g', 'MarkerSize', 8, 'MarkerEdgeColor','g'); 
 %subplot(1,2,2);
%  %surf (frFit);
%  pcolor(frFit); shading flat; colormap hot; title('1x 2D fit, psf fixed'); axis square ; hold on;
%  plot(x3+0.6,y3+0.6,'o','MarkerFaceColor', 'b', 'MarkerSize', 8, 'MarkerEdgeColor','b');
%  plot(x4+0.4,y4+0.4,'ro','MarkerFaceColor', 'g', 'MarkerSize', 8, 'MarkerEdgeColor','g'); 
%---------------------------------------------- 
 %[~]=ginput(1);
 %close(gcf);
end