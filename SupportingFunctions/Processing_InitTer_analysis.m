function Processing_InitTer_analysis(exp)
%Load database, analyze start and stop
%JacobKers 2012

close all
actions.savedata=1;         %default=1 (new analysis)
actions.loaddatabase=1;     %default=1 (new analysis)
actions.plot=0;             %default=0 (quick run)

if nargin<1, exp='TEST';end
initval=A001_Images_Set_Experiment(exp);
rnge=2*initval.extension;
outname=strcat(initval.basepath,initval.outname);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------
[~,chan_no]=size(S);

for ch=1:chan_no  %for each channel
chan_no-ch
[~,repno]=size(S(ch).channels.ReplicationCluster);
for lrp=1:repno  %for each bacterium     
%Get spot data---------------
time=S(ch).channels.ReplicationCluster(lrp).PosKyTracCom.frames_ext;
spots=S(ch).channels.ReplicationCluster(lrp).FluoPropsGen.fluospotscontent_ext;

tinitC=S(ch).channels.RepClicks(lrp).PosClick.firstframe;
tterC=S(ch).channels.RepClicks(lrp).PosClick.lastframe;

%select init and ter sections---------------
time_init=time(1:rnge);
spots_init=spots(1:rnge);

time_ter=time(end-rnge+1:end);
spots_ter=spots(end-rnge+1:end);

%padding to kill contribution of precessing and succesing replication spots
spots_initpad=spots_init;
[minval0,minidx0]=min(spots_initpad);
spots_initpad(1:minidx0)=minval0;  %separate init section

spots_terpad=spots_ter;
[minval1,minidx1]=min(spots_terpad); %separate ter section
spots_terpad(minidx1:end)=minval1;

spotspad=spots;
if 1
spotspad(1:minidx0)=minval0;
spotspad(end-rnge+minidx1:end)=minval1;
end

%Perform Stepfits--------------------------------------------------------
lsp=length(spotspad);
spl=Split2(spotspad,0,rnge); ixnitS=spl(2); 
tinitS=time(ixnitS);
spl=Split2(spotspad,lsp-rnge,lsp); ixterS=spl(2); 
tterS=time(ixterS);

%normalize time axis on results stepfit-------------------------
time_norm=(time-tinitS)/(tterS-tinitS);
time_init_norm=(time_init-tinitS)/(tterS-tinitS);
time_ter_norm=(time_ter-tinitS)/(tterS-tinitS);

%Store results---------------------------------------------------
S(ch).channels.ReplicationCluster(lrp).Cycle.normtime_ext=time_norm;
%This is the aligned time axis

S(ch).channels.ReplicationCluster(lrp).Cycle.Init.idx=ixnitS;
S(ch).channels.ReplicationCluster(lrp).Cycle.Init.time=tinitS;
S(ch).channels.ReplicationCluster(lrp).Cycle.Ter.idx=ixterS;
S(ch).channels.ReplicationCluster(lrp).Cycle.Ter.time=tterS;

S(ch).channels.ReplicationCluster(lrp).FluoPropsGen.fluospotscontentpadd_ext=spotspad;
%This is the fluorescent intensity with the contribution of former and
%later spots suppressed

if actions.plot
%plot menu------------------------------------
markY=[-0.2*max(spots) 1.2*max(spots)];
markXinit=([tinitC tinitC]-tinitS)/(tterS-tinitS);
markXter=([tterC tterC]-tinitS)/(tterS-tinitS);
markTinit=([tinitS tinitS]-tinitS)/(tterS-tinitS);
markTter=([tterS tterS]-tinitS)/(tterS-tinitS);


plot(time_norm,spots,'o-'); hold on;
plot(time_init_norm,spots_init,'go-');
plot(time_ter_norm,spots_ter,'go-')
plot(time_init_norm,spots_initpad,'ro-');
plot(time_ter_norm,spots_terpad,'ko-')

%Add cycle markers
plot(markXinit,markY, 'k-');
plot(markTinit,markY, 'r-');
plot(markXter,markY, 'k-');
plot(markTter,markY, 'r-');
pause(3);
close(gcf);
%[~]=ginput(1);
end
end
end
if actions.savedata
    save(outname,'S', '-append');
end
display('done');


function spl=Split2(rij,i1,i2)
	%this function adresses a one-dim array 'rij'in a specific segment
	%and determines the best step-fit there
	window=i2-i1;
	if window>2
        Chisq=(1:window-1)*0;
        for t=2:window-2;
            left=rij(i1+1:i1+t);
            right=rij(i1+t+1:i2);
            dcleft=mean(left);
            dcright=mean(right);
            Chisq(t)=(sum((left-dcleft).^2)+sum((right-dcright).^2))/(window-1);
        end
        Chisq(1)=Chisq(2)+1;
        Chisq(window-1)=Chisq(window-2)+1;
        [g,h]=min(Chisq);
        r=rij(i1+h+1:i2);
        l=rij(i1+1:i1+h);   
        stp=mean(r)-mean(l);           
        rank1=abs(stp)*sqrt(window);        %expected rel.step accuracy relative to background noise
        rank2=abs(stp)*sqrt(window/g);       %measured rel. step accuracy  
        spl=[sqrt(g),h+i1,stp,rank1, rank2];      %minimum Chi, index, stepsize,  step*srqt(N), step*sqrt(N/Var)
	else
        spl=[0,0,0,0,0];
	end