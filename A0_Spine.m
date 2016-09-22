clear
clc

%% Define experiment and user

exp ='Exp001_DnaN_TUS_dif_01092016_M';
user = 'Mark';

%% Prepare raw images (roicotrasca)

A010_Prepare_Images_After_NIS(user,exp);

%% Save Images to Workspace

A020_Images_SaveImageWorkSpace(user,exp);

%% Find Drift vector

A030_Images_FindDriftVector(user,exp)
close all

%% Kymomaker

initval = A040_Images_KymoMaker(user,exp);
WorkspaceOutName = initval.WorkspaceOutName;
close all

%% Click through cell cycles

A060_RepliCluster00_Click(user,exp,WorkspaceOutName)
close all

%% Processing

A100_Processing_ReplicationAutoShell(user,exp,WorkspaceOutName)

%% Save Bacpics to workspace

[Bacpics, MDchans, FLchans, Cells] = A110_Get_Bacpics(user,exp);
initval=A001_Images_Set_Experiment(user,exp);

for chan = 1:FLchans;
    IPTPvalue(chan) = LionfitcontrolUI_Kymo(initval,Bacpics{1,chan},Cells(chan),chan);
end


%% Lionfit

initval.Lionpath = strcat(initval.Kymopath,initval.OSslash,'LionFit',initval.OSslash);
addpath(genpath(strcat(initval.Lionpath,'150917V')));
addpath(genpath(strcat(initval.Lionpath,'gaussmlev2')));

GaussFitSimedit(user,exp)