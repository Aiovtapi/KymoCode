function Processing_Collect_DataBases(exp,secondexp)
%Load databases and save them in  common databases 'S' (for automatic
%processing results) and 'M'(anual) for user inputs (clicked positions
%etc.)

if nargin<1, exp='CM_DnaN_37Deg_Series1002';
end
%%------------------------------------
resaveclickdata=1;  %default=1 (new analysis)
%%%------------------------------------

initval=A001_Images_Set_Experiment(exp);


%[chans,~]=size(initval.nms) %think this was a mistake
[~,chans]=size(initval.nms)
chans

for i=1:chans
infi=strcat(initval.basepath,initval.nms{i});
buf=load(infi);
initval=A001_Images_Set_Experiment(exp);  %just to be sure
S(i).channels.initval=buf.initval;
S(i).channels.RepClicks=buf.RepClicks;

S(i).channels.kymo_FL=buf.kymo_FL;
S(i).channels.kymo_BF=buf.kymo_BF;
S(i).channels.chanstk_BF=buf.chanstk_BF;
S(i).channels.chanstk_FL=buf.chanstk_FL;
S(i).channels.ReplicationCluster=buf.ReplicationCluster;


M(i).channels.initval=buf.initval;
M(i).channels.endpoints=buf.endpoints;
M(i).channels.presets=buf.presets;
M(i).channels.RepClicks=buf.RepClicks;

[~,Nrep]=size(S);
for j=1:Nrep
 M(i).channels.RepClicks(j).accepted=1;
end

end

%-----------------------------------------------
%Saving. Note that the 'M' Database is NOT standard rewritten. This is
%because it contains manual input from various analysis stages (clicking
%bacterial cycles, accept-reject runs)
outnameS=strcat(initval.basepath,initval.outname);
save(outnameS, 'M','S');

outnameM=strcat(initval.basepath,initval.outname_usr);
if resaveclickdata, save(outnameM, 'M', '-append');end  
%only after re-clicking.watch out with this one!
disp('done');
%------------------------------------------------------
