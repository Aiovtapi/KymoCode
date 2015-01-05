function Processing_Collect_DataBases_2FLchan(exp,exp_FLchan2)
%Load databases and save them in  common databases 'S' (for automatic
%processing results) and 'M'(anual) for user inputs (clicked positions
%etc.)
close all;
if nargin<1, 
    exp='A_CM_DnaXDnaN_DualColour_Col002_DnaNSignal';
    exp_FLchan2='B_CM_DnaXDnaN_DualColour_Col002_DnaXSignal';
end
%%------------------------------------
resaveclickdata=1;  %default=1 (new analysis)
%%%------------------------------------

initval=A001_Images_Set_Experiment(exp);
initval_chan2=A001_Images_Set_Experiment(exp_FLchan2);

%[chans,~]=size(initval.nms) %think this was a mistake
[chans,~]=size(initval.nms)
chans

for i=1:chans
infi=strcat(initval.basepath,initval.nms{i});  %main channel (DnaN)
buf=load(infi);

infi2=strcat(initval.basepath,initval_chan2.nms{i});  %secondary channel (DnaX)
buf2=load(infi2);

initval=A001_Images_Set_Experiment(exp);  %just to be sure
initval_chan2=A001_Images_Set_Experiment(exp_FLchan2);

S(i).channels.initval=buf.initval;


S(i).channels.kymo_FL=buf.kymo_FL;
S(i).channels.kymo_FL2=buf2.kymo_FL;

S(i).channels.kymo_BF=buf.kymo_BF;
S(i).channels.chanstk_BF=buf.chanstk_BF;
S(i).channels.chanstk_FL=buf.chanstk_FL;
S(i).channels.chanstk_FL2=buf2.chanstk_FL;

if isfield(buf, 'ReplicationCluster')
    S(i).channels.ReplicationCluster=buf.ReplicationCluster;
end

S(i).channels=orderfields(S(i).channels);
[~,Nrep]=size(S);

M(i).channels.initval=buf.initval;
M(i).channels.endpoints=buf.endpoints;
M(i).channels.presets=buf.presets;

if isfield(buf, 'RepClicks')
    M(i).channels.RepClicks=buf.RepClicks;
    S(i).channels.RepClicks=buf.RepClicks;
    for j=1:Nrep
    M(i).channels.RepClicks(j).accepted=1;
    end
end


    if 1
    figure;
    subplot(1,3,1); pcolor(S(i).channels.kymo_FL); shading flat; colormap hot; title('DnaN');
    subplot(1,3,2); pcolor(S(i).channels.kymo_FL2); shading flat; colormap hot;title('DnaX');
    subplot(1,3,3); pcolor(S(i).channels.kymo_BF); shading flat; colormap hot;;title('Brightfield');
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
