function [pdf] = SingleDpdf(x,D)
%DOUBLEDPDF Summary of this function goes here
%   Detailed explanation goes here
n=4;
pdf=((n/D)^n*x^(n-1).*exp(-n*x/D))/(factorial(n-1));

end

