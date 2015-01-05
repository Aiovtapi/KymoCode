function E = TripleTwoDGaussexpected_Fix1PSF(x,y,psf,params)
% The expected counts per pixel.
[r,c]=size(x);

E = abs(params (6)).*twoDGauss_Symm_FixPSF (x,y,params (1),params (3),psf) +...  %peak 1
    abs(params (7)).*twoDGauss_Symm_FixPSF (x,y,params (2),params (4),psf) +...  %peak 2
    abs(params (5))*1/(2*pi*c/2*r/2)*exp(-((x-c/2).^2/(2*(c/2)^2)+(y-r/2).^2/(2*(c/2)^2))); %shallow background as elliptical shape
end

function p = twoDGauss_Symm_FixPSF(x,y,ux,uy,psf)
% Double 2D Gaussian; symmetric and fixed psf, normalized
p =   1/(2*pi*psf.^2)*exp (-( (x-ux).^2./(2*psf^2)+(y-uy).^2./(2*psf^2) ));
end