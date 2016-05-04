function Processing_Collect_DataBases(exp,ColourIdx,DnaNIdx)
%Load databases and save them in  common databases 'S' (for automatic
%processing results) and 'M'(anual) for user inputs (clicked positions
%etc.)

if nargin<1, exp='CM_DnaN_37Deg_Series1002';
end
%%------------------------------------
resaveclickdata=1;  %default=1 (new analysis)
%%%------------------------------------

initval=A001_Images_Set_Experiment(exp);

chans=initval.channelno;

for ch=1:chans
    DnaNpath=strcat(initval.basepath,initval.nms{ch}{DnaNIdx});
    bof=load(DnaNpath);
infi=strcat(initval.basepath,initval.nms{ch}{ColourIdx});
buf=load(infi);
initval=A001_Images_Set_Experiment(exp);  %just to be sure
S(ch).channels.initval=buf.initval;
S(ch).channels.RepClicks=bof.RepClicks;

S(ch).channels.kymo_FL=buf.kymo_FL;
S(ch).channels.kymo_BF=buf.kymo_BF;
S(ch).channels.chanstk_BF=buf.chanstk_BF;
S(ch).channels.chanstk_FL=buf.chanstk_FL;
S(ch).channels.ReplicationCluster=bof.ReplicationCluster;


M(ch).channels.initval=buf.initval;
M(ch).channels.endpoints=buf.endpoints;
M(ch).channels.presets=buf.presets;
M(ch).channels.RepClicks=bof.RepClicks;

[~,Nrep]=size(S);
for j=1:Nrep
 M(ch).channels.RepClicks(j).accepted=1;
end

end

%-----------------------------------------------
%Saving. Note that the 'M' Database is NOT standard rewritten. This is
%because it contains manual input from various analysis stages (clicking
%bacterial cycles, accept-reject runs)
% 
% if Doneclick;
%     outnameS=strcat(initval.basepath,initval.outname,'2');
% else
%     outnameS=strcat(initval.basepath,initval.outname);
% end
outnameS=strcat(initval.basepath,initval.outname{ColourIdx});

save(outnameS, 'M','S');

outnameM=strcat(initval.basepath,initval.outname_usr);
if resaveclickdata, save(outnameM, 'M', '-append');end  
%only after re-clicking.watch out with this one!
disp('done');
%------------------------------------------------------
