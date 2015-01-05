function p = twoDGauss (x,y,ux,uy,s1,s2)
% Double 2D Gaussian
p = 1/(2.*pi.*s1.*s2) .* exp (-( (x-ux).^2./(2.*s1.^2)+(y-uy).^2./(2.*s2.^2) ));
end