function Sout=DataBase_Reorganizer(S);

%S.channels.Division:-------------------------

%           name: 2048
%         frames: [1x28 double]
%           left: [1x28 double]
%          right: [1x28 double]
%           fate: 'divided'
%        firstfr: 286
%         lastfr: 313
%      firstleft: 14
%     firstright: 34
%       lastleft: 17
%      lastright: 53
%      linkedrep: 24

%S.channels.ReplicationCluster------------------:
%           name: 2048
%         frames: [1x23 double]
%       clickpos: [1x23 double]
%       trackpos: [23x1 double]
%           fate: 'disassembled'
%       spot1pos: [23x1 double]
%       spot2pos: [23x1 double]
%       firstpos: 24
%        lastpos: 32
%     firstframe: 280
%      lastframe: 302
%       spot1amp: [23x1 double]
%       spot1sig: [23x1 double]
%       spot2amp: [23x1 double]
%       spot2sig: [23x1 double]
%      linkedbac: 27
%--------------------------------


[~,chans]=size(S);
clear Sout
for ch=1:chans  %for each channel
        Sout(ch).channels.initval=S(ch).channels.initval;
        Sout(ch).channels.kymo_FL=S(ch).channels.kymo_FL;
        Sout(ch).channels.kymo_BF=S(ch).channels.kymo_BF;
        Sout(ch).channels.chanstk_BF= S(ch).channels.chanstk_BF;
        Sout(ch).channels.chanstk_FL=S(ch).channels.chanstk_FL;
  
    [bac,~]=size(S(ch).channels.Division);
    
    for bc=1:bac
        bc
         %basic labels
         Sout(ch).channels.Division(bc).name=S(ch).channels.Division(bc).name;
         Sout(ch).channels.Division(bc).fate=S(ch).channels.Division(bc).fate;
         Sout(ch).channels.Division(bc).linkedrep=S(ch).channels.Division(bc).linkedrep;
         
         %positions: cycle limits, manual points
         Sout(ch).channels.Division(bc).PosClick.firstfr=S(ch).channels.Division(bc).firstfr(:);  
         Sout(ch).channels.Division(bc).PosClick.lastfr=S(ch).channels.Division(bc).lastfr(:);
         Sout(ch).channels.Division(bc).PosClick.firstleft=S(ch).channels.Division(bc).firstleft(:);
         Sout(ch).channels.Division(bc).PosClick.firstright=S(ch).channels.Division(bc).firstright(:);
         Sout(ch).channels.Division(bc).PosClick.lastleft=S(ch).channels.Division(bc).lastleft(:);
         Sout(ch).channels.Division(bc).PosClick.lastright=S(ch).channels.Division(bc).lastright(:);
         
         %positions: tracked from first kymograph-based analysis
         Sout(ch).channels.Division(bc).PosKyTrac.frames=S(ch).channels.Division(bc).frames;
         Sout(ch).channels.Division(bc).PosKyTrac.left=S(ch).channels.Division(bc).left;
         Sout(ch).channels.Division(bc).PosKyTrac.right=S(ch).channels.Division(bc).right;
         
    end
    [rep,~]=size(S(ch).channels.ReplicationCluster);
    for rp=1:rep
        rp;
         %basic labels
         Sout(ch).channels.ReplicationCluster(rp).name=S(ch).channels.ReplicationCluster(rp).name;
         Sout(ch).channels.ReplicationCluster(rp).fate=S(ch).channels.ReplicationCluster(rp).fate;
         Sout(ch).channels.ReplicationCluster(rp).linkedbac=S(ch).channels.ReplicationCluster(rp).linkedbac;
         
         %positions: cycle limits, manual points
         Sout(ch).channels.ReplicationCluster(rp).PosClick.firstframe=S(ch).channels.ReplicationCluster(rp).firstframe;
         Sout(ch).channels.ReplicationCluster(rp).PosClick.lastframe=S(ch).channels.ReplicationCluster(rp).lastframe;
         Sout(ch).channels.ReplicationCluster(rp).PosClick.firstpos=S(ch).channels.ReplicationCluster(rp).firstpos;
         Sout(ch).channels.ReplicationCluster(rp).PosClick.lastpos=S(ch).channels.ReplicationCluster(rp).lastpos;      
         
          %positions: tracked from first kymograph-based analysis: compos
          Sout(ch).channels.ReplicationCluster(rp).PosKyTracCom.frames=S(ch).channels.ReplicationCluster(rp).frames;
          Sout(ch).channels.ReplicationCluster(rp).PosKyTracCom.clickpos=S(ch).channels.ReplicationCluster(rp).clickpos;
          Sout(ch).channels.ReplicationCluster(rp).PosKyTracCom.trackpos=S(ch).channels.ReplicationCluster(rp).trackpos;

          %positions: tracked from first kymograph-based analysis: double gaussian
           Sout(ch).channels.ReplicationCluster(rp).PosKyTracGauss.spot1pos=S(ch).channels.ReplicationCluster(rp).spot1pos;
           Sout(ch).channels.ReplicationCluster(rp).PosKyTracGauss.spot1amp=S(ch).channels.ReplicationCluster(rp).spot1amp;
           Sout(ch).channels.ReplicationCluster(rp).PosKyTracGauss.spot1sig=S(ch).channels.ReplicationCluster(rp).spot1sig;
           Sout(ch).channels.ReplicationCluster(rp).PosKyTracGauss.spot2pos=S(ch).channels.ReplicationCluster(rp).spot2pos;
           Sout(ch).channels.ReplicationCluster(rp).PosKyTracGauss.spot2amp=S(ch).channels.ReplicationCluster(rp).spot2amp;
           Sout(ch).channels.ReplicationCluster(rp).PosKyTracGauss.spot2sig=S(ch).channels.ReplicationCluster(rp).spot2sig;
%--------------------------------              
    end
end
