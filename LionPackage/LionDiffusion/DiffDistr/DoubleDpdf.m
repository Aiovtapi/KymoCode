function [pdf] = DoubleDpdf(x,D1,D2)
%DOUBLEDPDF Summary of this function goes here
%   Detailed explanation goes here
n=4;
A=0.5;
pdf=A*((n/D1)^n*x.^(n-1).*exp(-n.*x./D1))/(factorial(n-1))+(1-A)*((n/D2)^n*x.^(n-1).*exp(-n.*x./D2))/(factorial(n-1));
end

