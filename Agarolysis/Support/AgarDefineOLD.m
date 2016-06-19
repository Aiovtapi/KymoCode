function init = AgarDefine(user)


init.difchan = 'CFP';

AgarUI(init)

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
%% Add paths

init.Agarpath = strcat(init.kymopath,'Agarolysis',init.OSslash);

addpath(init.Agarpath);
switch init.OSslash
    case '\'
        addpath(genpath(strcat(init.Agarpath,'oufti_windows')));
    case '/'
        addpath(genpath(strcat(init.Agarpath,'Oufti_source_code')));
end
addpath(strcat(init.kymopath,'LionFit',init.OSslash,'150917V'));
addpath(strcat(init.kymopath,'Agarolysis',init.OSslash,'Support'));


%%
init.bfimgname = 'BF.tif';

defaultset = questdlg('Load default dataset?','Menu','Yes','No','Yes');
switch defaultset
    case 'Yes'
        init.pcimgname = 'PC.tif';
        init.CFPimgname = '457-100ms-10mWo-300G.tif';
        init.YFPimgname = '515-100ms-50mWo-300G.tif';
        init.RFPimgname = '561-100ms-33mWo-300G.tif';
        init.CFPbeampath = strcat(init.kymopath,'BeamShape457.tif');
        init.YFPbeampath = strcat(init.kymopath,'BeamShape515.tif');
        init.RFPbeampath = strcat(init.kymopath,'BeamShape561.tif');
        init.meshfile = 'PC.mat';
        init.meshpath = strcat(init.datapath,init.meshfile);
        init.pcresize = 0.421;
        init.pctrans = [0,0];
        init.flresize = 1; 
        init.fltrans = [2,-63];
    case 'No'
        init.datapath = uigetdir('','Select data folder');
        init.datapath = strcat(init.datapath,init.OSslash);
        cd(init.datapath)
        init.pcimgname = uigetfile('*.tif','Select PC image');
        init.CFPimgname = uigetfile('*.tif','Select CFP image');
        init.YFPimgname = uigetfile('*.tif','Select YFP image');
        init.RFPimgname = uigetfile('*.tif','Select RFP image');
        [meshname,meshpath,~] = uigetfile('*.mat','Select oufti output for dataset');
        init.meshpath = strcat(meshpath,meshname);
        
        defaultbeam = questdlg('Load default beamshapes?','Menu','Yes','No','Yes');
        switch defaultbeam
            case 'Yes'
                init.CFPbeampath = strcat(init.kymopath,'BeamShape457.tif');
                init.YFPbeampath = strcat(init.kymopath,'BeamShape515.tif');
                init.RFPbeampath = strcat(init.kymopath,'BeamShape561.tif');
            case 'No'
                [a,b,~] = uigetfile('*.tif','Select CFP beamshape');
                init.CFPbeampath = strcat(b,a);
                [a,b,~] = uigetfile('*.tif','Select YFP beamshape');
                init.YFPbeampath = strcat(b,a);
                [a,b,~] = uigetfile('*.tif','Select RFP beamshape');
                init.RFPbeampath = strcat(b,a);
        end
        
        defaulttrans = questdlg('Use default translations?','Menu','Yes','No','Yes');
        switch defaulttrans
            case 'Yes'
                init.pcresize = 0.421;      % scaling factor for phase contrast 
                init.pctrans = [0,0];       % translation for phase contrast [x,y]
                init.flresize = 1;          % scaling factor for fluorescence
                init.fltrans = [2,-63];     % translation of fluorescence [x,y]
            case 'No'
                [ans1, ans2, ans3] = inputdlg({'PC scale factor','PC x translation','PC y translation'}...
                    ,'PC',1,{'0.421','0','0'},'on');
                [ans4, ans5, ans6] = inputdlg({'FL scale factor','FL x translation','FL y translation'}...
                    ,'FL',1,{'1','2','-63'},'on');
                init.pcresize = str2double(ans1);
                init.pctrans = [str2double(ans2),str2double(ans3)];
                init.flresize = str2double(ans4);
                init.fltrans = [str2double(ans5),str2double(ans6)];
        end
                
end
        
init.maxfile = 421;         % max amount of pictures allowed in stack
init.lioncropindex = 0;     % whether bacpics are cropped in lionfit
init.Extrabound = 4;        % extra boudaries added to bacpics
init.strelval = 8;          % disk radius for imdilate of bacpic mask
init.IPTP = 1;              % Intensity Peak threshold parameter
init.bacpath = strcat(init.datapath,'Bacpics',init.OSslash);

%% Selection of channels

options = {'CFP','YFP','RFP'};

viewchannels = listdlg('PromptString','Select channels to view:','SelectionMode','multiple',...
                'ListString',options);
            
for idx = 1:numel(viewchannels);
    switch viewchannels(idx)
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