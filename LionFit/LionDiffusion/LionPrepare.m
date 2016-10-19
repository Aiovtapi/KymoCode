function [Traj] = LionPrepare(N,LB,UB)
%LIONPREPARE prepare utrack trajectories for VB3 analysis

% N : numbers of experiments considered for analysis

% loadutrack is function that loads trajectories in variable 'tracksFinal'
% 'tracksFinal' is cell array, each index different utrack analysis.
% each cell contains output structure of utrack

% Measurements are usually separated by numbers, e.g. 1, 2, 3 .. 
% utrackResults are saved in these folder separately

% first copy these folders to main utrack folder in the results dir

curFolder=pwd;
prevFolder=fileparts(pwd);

for i=N
    if exist(strcat(curFolder,'/utrack/utrackResults',num2str(i)))==0
    utrackfolder=strcat(prevFolder,'/',num2str(i),'/utrackResults');
    copyfile(utrackfolder,strcat(curFolder,'/utrack/utrackResults',num2str(i)))
    end
end

% Then create cell with all trajectories concatentated.

tracksFinal=loadutrack(N);

% Take all trajectories and store them in VB3 format
% boundary conditions are applied for trajectory min (LB) max lengths (UB)

dim=2; %dimensionality

%Main trajectory cell:
F=utracktoVB3(tracksFinal,LB,UB,dim);

%concatenate trajectories from cell
sizeF=size(F,2);
Traj=cell(1,1);
for i=1:size(F,2)
    Traj=[Traj F{i}.Trajectories];
end

Traj=Traj(~cellfun('isempty',Traj));

if exist('Data/Traj.mat')==2
    delete('Data/Traj.mat')
    mkdir('Data')
    save('Data/Traj.mat','Traj');
elseif exist('Data/Traj.mat')==0
    mkdir('Data');
    save('Data/Traj.mat','Traj');
end

end

