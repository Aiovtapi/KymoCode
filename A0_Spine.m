clear all 
clc

%% Define experiment and user

exp ='Exp001_DnaN_TUS_dif_01092016_M';
user = 'MarkPC';

%% Prepare raw images (roicotrasca)

A010_Prepare_Images_After_NIS(user,exp);

%% Save Images to Workspace

A020_Images_SaveImageWorkSpace(user,exp);

%% Find Drift vector

A030_Images_FindDriftVector(user,exp)
close all

%% Kymomaker

A040_Images_KymoMaker(user,exp)

%% Click through cell cycles

A060_RepliCluster00_Click(user,exp)

%% Processing

A100_Processing_RepclclicationAutoShell(user,exp)