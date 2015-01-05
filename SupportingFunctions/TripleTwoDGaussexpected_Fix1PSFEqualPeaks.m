function E = TripleTwoDGaussexpected_Fix1PSFEqualPeaks(x,y,psf,params)
% The expected counts per pixel.
[r,c]=size(x);
a=c/2;b=r/2;

% E = params (6).*twoDGauss_Symm_FixPSF (x,y,params (1),params (3),psf) +...
%     params (7).*twoDGauss_Symm_FixPSF (x,y,params (2),params (4),psf) +...
%     params (5).*exp(-(((x-a)/a).^2+((y-b)/b).^2)); %independent peaks


E = params (6).*twoDGauss_Symm_FixPSF (x,y,params (1),params (3),psf) +...
    params (6).*twoDGauss_Symm_FixPSF (x,y,params (2),params (4),psf) +...
    params (5).*exp(-(((x-a)/a).^2+((y-b)/b).^2)); %equal peaks
end