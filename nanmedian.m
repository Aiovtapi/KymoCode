function [Y]=nanmedian(X) 
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
X(isnan(X))=0;
X=nonzeros(X);
Y=median(X);
end