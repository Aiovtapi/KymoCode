%Processing_FetchDataFromDataBase
%Load database, 
%JacobKers 2012----------------------------------
close all
%exp='TEST'; %See 'Images_Set_Experiment(exp)';
exp='001_DnaN_TUS_dif_30122014_DnaNsignal_test';
minsperframe=2.5;

initval=A001_Images_Set_Experiment(exp);
outname=strcat(initval.basepath,initval.outname);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
load(outname,'S');
load(outname_usr,'M');
%------------------------------------------------------------------

[~,chan_no]=size(S);

for ch=1:chan_no  %for each channel
chan_no-ch   
Div=S(ch).channels.AutoDivision;
Fluo1=S(ch).channels.ReplicationCluster;
% Fluo2=S(ch).channels.SecondFluoCluster;
[~,bacno]=size(Div);
for bc=1:bacno  %for each bacterium 
 ok=S(ch).channels.AutoDivision(bc).accepted  %manual accept/reject
if ok
    %Get fluorescent values-------------------------------------
    frames=minsperframe*S(ch).channels.ReplicationCluster(bc).PosKyTracCom.frames;  %0) Select the time axis (in frames)
    switch 1
        case 1
        FL_A=S(ch).channels.ReplicationCluster(bc).Pos2DPreTrac.contentallspots'; %1) Select total count
        FL_B=S(ch).channels.SecondFluoCluster(bc).Pos2DPreTrac.contentallspots'; %2) select  spots count
         case 2
        FL_A=S(ch).channels.ReplicationCluster(bc).Pos2DFinTrac.contentallspots'; %1) Select total count
        FL_B=S(ch).channels.SecondFluoCluster(bc).Pos2DFinTrac.contentallspots'; %2) select  spots count
        case 3
        FL_A=S(ch).channels.ReplicationCluster(bc).FluoPropsFin.content_spots';
        FL_B=S(ch).channels.SecondFluoCluster(bc).FluoPropsFin.content_spots';
    end
    %just plot it
    figure;
    subplot(2,1,1); plot(frames,FL_A,'o-');
    title(strcat('Chan',num2str(ch), 'BacName',num2str(bc)));
    legend('spots DnaN');
    subplot(2,1,2);plot(frames,FL_B, 'r-o');
    legend('spotsDnaX');
    
    
    
    [~]=ginput(1);
    close(gcf);
end
end
end
