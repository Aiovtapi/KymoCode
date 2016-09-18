function [Y] = nanmin(X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
X(isnan(X))=0;
X=nonzeros(X);
Y=min(X);
end

