function initval=Exp001_DnaN_TUS_dif_30122014_M(initval, user) 

initval.viewchannel = 'RFP';

switch user
      case 'Jacob', 
%           initval.basepath='D:\jkerssemakers\My Documents\BN_ND_Data_Recent Master\BactRepl_Charl\20132511_NewExps\GammaMinus\Series1\';
%           initval.BFdatapath='Brightfield_AVGSub\';
%           initval.FiguresFolder='Figures\';
%           initval.FLdatapath='Fluorescence_Aligned\';
      case 'Charl_Local' ,
          %initval.basepath='/Volumes/LittleMonster/PhD_Docs/Gamma_minus_Experiments/Data/20130724-DnaN-Dif-Gamma-ve/Tiffs/Series1/';
%           initval.basepath='/Users/mcmoolman/Documents/Publications/BetaClampStoichiometry_and_Dynamics/20140117_Data/DnaX_YPet_AND_mCherry_DnaN/20130621_DnaX_DnaN_Colocalization/DnaX_mYPet_mCherry_DnaN_col002/';
% 
%           initval.BFdatapath='Brightfield/';
%           initval.FiguresFolder='Figures/';
%           initval.FLdatapath='Fluorescence/DnaX_YPet_Fluorescence_Aligned/';
    case 'Roy'
        initval.basepath='/Users/rleeuw/Work/Data/141230_dnaN_dif_tus/dnaN_dif_tus/';
        initval.BFdatapath='Brightfield/';
        initval.FiguresFolder='Figures/';
        initval.FLpath='Figures/';
    case 'RoyPC'
%         initval.basepath='D:\rleeuw\Data\141230_dnaN_dif_tus\dnaN_dif_tus_40msExpTime_5minAcqTimeYFP_30msExpCFP_002_C1\';
%         initval.BFdatapath='Brightfield\';
%         initval.FiguresFolder='Figures\';
%         initval.FLdatapath='Fluorescence\RFP\';

    case 'MarkPC'
        initval.basepath='D:\Users\water_000\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\';
        initval.BFdatapath='Brightfield\';
        initval.FiguresFolder='Figures\';
        initval.FLpath='Fluorescence\';
        
    case 'Mark'
        initval.basepath='D:\Users\water\OneDrive\Documents\BEP\Data\141230_dnaN_dif_tus\';
        initval.BFdatapath='Brightfield\';
        initval.FiguresFolder='Figures\';
        initval.FLpath='Fluorescence\';
        
end

initval.FLdatapath=[initval.FLpath,initval.viewchannel,'\'];

initval.BFfiletemplate='BF1.tif';
initval.FLfiletemplate=[initval.viewchannel,'_R_I_Corrected.tif'];
initval.maxfile=421;

%used for kymograph processing-------------------------
initval.channeldistance=32.0; %in pixels
initval.channelno=1; %number of channels
initval.perpadjust=0; %degrees
initval.kymoangle=47;  %CCW
initval.kymolength=150;  %in pixels
initval.kymohwidth=10;  %in pixels
initval.correctdrift=1;
initval.entranceoffset=0;  %12
%----------------------------------------

%Some common settings
initval.driftfile='driftvector.txt'; %used for drift correction
initval.sumHW=3;  %halfwidth of summaton area, w=2*sumHW+1
initval.spotRoiHW=9;  %This is a halfwidth used for fluorescence analysis 
initval.psf=1.3;    %point spead function 
initval.spotregionhalfwidth=9;  %This sets the area for detailed spot analysis
initval.startzoom=20;
initval.skip2Danalysis=1;


%Used for database collection----------------
%initval.nms={'ExpCM_DnaXDnaN_DualColour_Col002_DnaXSignalChan_x209' ;'ExpCM_DnaXDnaN_DualColour_Col002_DnaXSignalChan_x231'; 'ExpCM_DnaXDnaN_DualColour_Col002_DnaXSignalChan_x273'; 'ExpCM_DnaXDnaN_DualColour_Col002_DnaXSignalChan_x295'; 'ExpCM_DnaXDnaN_DualColour_Col002_DnaXSignalChan_x317'; 'ExpCM_DnaXDnaN_DualColour_Col002_DnaXSignalChan_x339'} ; 
%initval.nms={'ExpCM_Dn0aN-Dif-Gamma-ve-Position1_Series1Chan_x126' ;'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x146'; 'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x173'; 'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x219'; 'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x245'; 'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x291'; 'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x314'; 'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x382'; 'ExpCM_DnaN-Dif-Gamma-ve-Position1_Series1Chan_x405'} ; 

initval.nms={'Exp001_DnaN_TUS_dif_30122014_DnaNsignal_MChan_x228'};


initval.outname='DnaN_TUS_dif';
initval.outname_usr='DnaN_TUS_dif_UserInputs';
%Used for first analysis of fluorescence---------
initval.extension=7;
initval.BW_Threshold=1;

%used for clicking-------------------------------------
initval.zoom=100;
initval.estimateddoublingtime=48; %in frames   
%initval.zoom=1.3*initval.estimateddoublingtime;



