function [init] = AgarDefine(user,init)

if nargin < 2
    init.viewchan = 'CFP';
end

init.difchan = 'CFP';

switch user
    case 'Mark'
        init.OSslash = '\';
        init.kymopath = 'C:\Users\water\Documents\GitHub\KymoCode\';
        init.datapath = 'C:\Users\water\Documents\GitHub\Data\DnaN_dif_Tus_AgarPad\';
    case 'MarkPC'
        init.OSslash = '\';
        init.kymopath = 'D:\Users\water_000\Documents\GitHub\KymoCode\';
        init.datapath = 'D:\Users\water_000\Documents\GitHub\Data\DnaN_dif_Tus_AgarPad\';
end

init.bfimgname = 'BF.tif';
init.pcimgname = 'PC.tif';
init.CFPimgname = 'Simul_457-100ms-10mWo-300G.tif';
init.YFPimgname = '515-100ms-10mWo-300G.tif';
init.RFPimgname = '561-100ms-10mWo-300G.tif';
init.CFPbeamshape = 'BeamShape457.tif';
init.YFPbeamshape = 'BeamShape515.tif';
init.RPFbeamshape = 'BeamShape561.tif';

init.meshesfile = 'Simul_PC.mat';

init.maxfile = 421;

init.pcresize = 0.421;
init.pctrans = [0,0];
init.flresize = 1;
init.fltrans = [2,-63];

init.lioncropindex = 0;

%% Non edits from now

init.Agarpath = strcat(init.kymopath,'Agarolysis',init.OSslash);
init.bacpath = strcat(init.datapath,'Bacpics',init.OSslash);

addpath(init.Agarpath);
switch init.OSslash
    case '\'
        addpath(genpath(strcat(init.Agarpath,'oufti_windows')));
    case '/'
        addpath(genpath(strcat(init.Agarpath,'Oufti_source_code')));
end
addpath(strcat(init.kymopath,'LionFit',init.OSslash,'150917V'));

switch init.viewchan
    case 'CFP'
        init.flimgname = init.CFPimgname;
        init.beamshape = init.CFPbeamshape;
    case 'YFP'
        init.flimgname = init.YFPimgname;
        init.beamshape = init.YFPbeamshape;
    case 'RFP'
        init.flimgname = init.RFPimgname;
        init.beamshape = init.RFPbeamshape;
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