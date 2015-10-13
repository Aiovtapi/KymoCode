function [Y] = nanmean(X)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
X(isnan(X))=0;
X=nonzeros(X);
Y=mean(X);
end

