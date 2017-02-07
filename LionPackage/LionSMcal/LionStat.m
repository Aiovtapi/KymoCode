clc
clear all

%% Load Results

TraceNumbers=[9 11];

init.OSslash = '/';

fprintf('Select Final Results Folder');
init.resultspath = uigetdir(pwd,'Select Results Folder');
    
if init.resultspath == 0;
    init.resultspath = '/Users/rleeuw/Work/Data/170111_Tus-SMcal/gain100/FinalResults';
end

init.resultspath = strcat(init.resultspath,init.OSslash);

Results=cell(length(TraceNumbers),1);

w=1;
for i=TraceNumbers
    Results{w}=load(strcat(init.resultspath,'SMCResult',num2str(i)));
    w=w+1;
end

%% Collect all integrated intensity data
IItot=[];

for i=1:size(Results,1)
    IItot=[IItot;Results{i}.SMCResult.TraceVal];
end

%% plot histogram

hist(IItot,50)

