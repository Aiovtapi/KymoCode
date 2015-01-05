function p = twoDGauss_Symm_FixPSF(x,y,ux,uy,psf)
% Double 2D Gaussian; symmetric and fixed psf, normalized
p =   1/(2*pi*psf.^2)*exp (-( (x-ux).^2./(2.*psf.^2)+(y-uy).^2./(2.*psf.^2) ));
end