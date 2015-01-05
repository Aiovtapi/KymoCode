function E = DoubleTwoDGauss_Fix1PSF(x,y,psf,params)
% The expected counts per pixel.
[r,c]=size(x);

E = abs(params (5)).*twoDGauss_Symm_FixPSF (x,y,params (1),params (3),psf) +...  %peak 1
    abs(params (6)).*twoDGauss_Symm_FixPSF (x,y,params (2),params (4),psf);  %peak 2
end

function p = twoDGauss_Symm_FixPSF(x,y,ux,uy,psf)
% Double 2D Gaussian; symmetric and fixed psf, normalized
p =   1/(2*pi*psf.^2)*exp (-( (x-ux).^2./(2*psf^2)+(y-uy).^2./(2*psf^2) ));
end