function init = AgarDefine

%% Other parameters

init.maxfile = 421;         % max amount of pictures allowed in stack
init.lioncropindex = 0;     % whether bacpics are cropped in lionfit
init.Extrabound = 4;        % extra boudaries added to bacpics
init.strelval = 4;          % disk radius for imdilate of bacpic mask
init.difchan = 'CFP';
init.channels = {'CFP', 'YFP', 'RFP'}; % maximum of 3
init.IPTP = 1.5;  % Intensity Peak threshold parameter

%% Load UI

[init, ddata, dtrans, dbeam] = AgarUI(init);

%% Add paths

addpath(init.Agarpath);
switch init.OSslash
    case '\'
        addpath(genpath(strcat(init.Agarpath,'oufti_windows')));
    case '/'
        addpath(genpath(strcat(init.Agarpath,'Oufti_source_code')));
end
addpath(strcat(init.kymopath,'LionFit',init.OSslash,'150917V'));
addpath(genpath(strcat(init.kymopath,'LionFit',init.OSslash,'gaussmlev2')));
addpath(strcat(init.kymopath,'Agarolysis',init.OSslash,'Support'));

%% Default values

if ddata == 1
    init.datapath = 'D:\Users\water_000\Documents\GitHub\Data\Target Data\Agar Data\1\';
    init.PCimgname = 'PC.tif';

    init.CFPimgname = '457-100ms-10mWo.tif';
    init.YFPimgname = '515-100ms-50mWo.tif';
    init.RFPimgname = '561-100ms-33mWo.tif';
    init.meshfile = 'PC1.mat';
    init.meshpath = strcat(init.datapath,init.meshfile);

    init.CFPimgname = '457-100ms-10mWo-300G.tif';
    init.YFPimgname = '515-100ms-50mWo-300G.tif';
    init.RFPimgname = '561-100ms-33mWo-300G.tif';
    init.meshfile = 'PC.mat';
    init.meshpath = strcat(init.datapath,init.OSslash,init.meshfile);
end

init.bacpath = strcat(init.datapath,'Bacpics',init.OSslash);

if dtrans == 1;
    init.pcresize = 0.421;
    init.pctrans = [2,-2];
    init.flresize = 1; 
    init.fltrans = [2,-63];
end
    
if dbeam == 1;
    init.CFPbeampath = strcat(init.kymopath,'BeamShape457.tif');
    init.YFPbeampath = strcat(init.kymopath,'BeamShape515.tif');
    init.RFPbeampath = strcat(init.kymopath,'BeamShape561.tif');
end
       
%% Selection of channels
            
for idx = 1:numel(init.viewchannels);
    switch init.viewchannels(idx)
        case 1
            init.flimgname{idx} = init.CFPimgname;
            init.beampath{idx} = init.CFPbeampath;
        case 2
            init.flimgname{idx} = init.YFPimgname;
            init.beampath{idx} = init.YFPbeampath;
        case 3
            init.flimgname{idx} = init.RFPimgname;
            init.beampath{idx} = init.RFPbeampath;
    end
end

switch init.difchan
    case 'CFP'
        init.difimgname = init.CFPimgname;
    case 'YFP'
        init.difimgname = init.YFPimgname;
    case 'RFP'
        init.difimgname = init.RFPimgname;
end
end