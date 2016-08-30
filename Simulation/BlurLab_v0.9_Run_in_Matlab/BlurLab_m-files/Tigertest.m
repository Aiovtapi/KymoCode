clc

ini.SigmaR = 0.03;  % Variance for normal dis. used for change in R (displacement)
ini.SigmaA = 0.03;   % Variance for normal dis. used for change in A (angle)
ini.SigmaI = 0.03;   % Variance for normal dis. used for change in I (Intensity)

ini.MergeD = 3;   % Distance between two spots needed for merging
ini.CMerge = 0.2;   % Chance for a merging event (if within distance)

ini.CSplit = 0.02;   % Chance for a splitting event
ini.DSplit = 0.5;   % Change in displacement for one particle

ini.CDis = 0.01;       % Chance for a disappearence event

ini.Apois = 0.5;     % Number of spots appearences (poisson distribution)

ini.Nblink = 5;   % Percentage of blinking spots
ini.Tblink = 5;   % Average blinking duration (in frames)
ini.Cblink1 = 0.03; % Chance of becoming a blinking spot (new spots)
ini.Cblink2 = 0.3; % Chance of becoming a stable spot (blinking spots)
ini.MinI = 25;        % Minimal value of intensity

ini.tog_merging = 1;
ini.tog_split = 1;
ini.tog_spotdis = 1;
ini.tog_spotapp = 1;
ini.tog_blinking = 1;
ini.tog_dimer = 1;
ini.Tigfolder = uigetdir(pwd);

nframes = 50;
pts = 50;
meanI = 100;
D = 0.2;
Lx = 100;
Ly = 100;
Lz = 1;


Tigercreate(nframes,pts,meanI,D,Lx,Ly,Lz,ini)