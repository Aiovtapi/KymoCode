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
    case 'Roy'
        init.OSslash = '/';
        init.kymopath = '/Users/rleeuw/Work/DataAnalysis/201511_TusdifDnaN_Montage/';
        init.datapath = '/Users/rleeuw/Work/Data/160205_BN2384_and_Beam_Profiles/5/';
end

init.bfimgname = 'BF.tif';
% init.pcimgname = 'PC.tif';
% init.CFPimgname = '457-100ms-10mWo-300G.tif';
% init.YFPimgname = '515-100ms-50mWo-300G.tif';
% init.RFPimgname = '561-100ms-33mWo-300G.tif';
init.pcimgname = 'StackedPC.tif';
init.CFPimgname = 'Stacked457.tif';
init.YFPimgname = 'Stacked515.tif';
init.RFPimgname = 'Stacked561.tif';
init.CFPbeamshape = 'BeamShape457.tif';
init.YFPbeamshape = 'BeamShape515.tif';
init.RFPbeamshape = 'BeamShape561.tif';

init.meshesfile = 'StackedPC.mat';

init.maxfile = 421;

init.pcresize = 0.421;      % scaling factor for phase contrast 
init.pctrans = [0,0];       % translation for phase contrast [x,y]
init.flresize = 1;          % scaling factor for fluorescence
init.fltrans = [2,-63];     % translation of fluorescence [x,y]

init.lioncropindex = 0;     % whether bacpics are cropped in lionfit

init.Extrabound = 4;
init.strelval = 8;     % disk radius for imdilate of bacpic mask

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
addpath(strcat(init.kymopath,'Agarolysis',init.OSslash,'Support'));

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