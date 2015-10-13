function [Y]=nanstd(X)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
X(isnan(X))=0;
X=nonzeros(X);
Y=std(X);
end

