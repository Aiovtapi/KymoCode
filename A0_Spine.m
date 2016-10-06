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

if exist('WorkspaceOutName','var')
    A060_RepliCluster00_Click(user,exp,WorkspaceOutName)
else
    A060_RepliCluster00_Click(user,exp)
end
close all

%% Processing, producing bacpics

if exist('WorkspaceOutName','var')
    A100_Processing_ReplicationAutoShell(user,exp,WorkspaceOutName)
else
    A100_Processing_ReplicationAutoShell(user,exp)
end

%% Save Bacpics to workspace

[Bacpics, MDchans, FLchans, Cells] = Get_Bacpics(user,exp);


%%

initval=A001_Images_Set_Experiment(user,exp);
finddif = strcmp(initval.difchan,initval.viewchan);
viewFLchans = {initval.viewchan{finddif},initval.viewchan{~finddif}};
IPTP = 1;

for chan = 1:FLchans;
    IPTPvalue(chan) = LionfitcontrolUI(1,Bacpics{1,chan},viewFLchans,Cells(1),chan,'Kymo');
end


%% Lionfit

initval.Lionpath = strcat(initval.Kymopath,initval.OSslash,'LionFit',initval.OSslash);
addpath(genpath(strcat(initval.Lionpath,'150917V')));
addpath(genpath(strcat(initval.Lionpath,'gaussmlev2')));

if exist('IPTPvalue','var')
    GaussFitSimedit(user,exp,IPTPvalue)
else
    GaussFitSimedit(user,exp)
end
