function [A] = loadutrack(N)
%LOADUTRACK Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
    N=[];
end

mainstring='./utrack/utrackResults';
followup='/TrackingPackage/tracks/Channel_1_tracking_result.mat';

for i=1:length(N)
    
    Nstr=num2str(N(i));
    tracksFinal{i}=load(strcat(mainstring,Nstr,followup));
    
end

A=tracksFinal;
end

