%function A102_Processing_Find_Division_Times(exp)
function Processing_Find_Division_Times(exp,ColourIdx)
%Load database, 
%1) For every replicationcluster....
    %find the two daugthers
    %(if found) create 'left'and 'right' forks that start at termination and extend to the new
    %daugther init points. The new division should lie in the region encompassed by this fork and the subsequent region.
    %this division marks the end of this bacterium.
%We use a binary brightfield-kymograph to do this

%JacobKers 2012

close all
actions.savedata=1;         %default=1 (new analysis)
actions.loaddatabase=1;     %default=1 (new analysis)
actions.plot=1;             %default=0 (quick run)

if nargin<1, 
    exp='006_oriZ_dif_110515_DnaNsignal';
    %exp='CM_DnaN-Dif-Gamma-ve-Position1_Series1';
    %exp='CM_HFR1_JK' ; 
    %exp='2a';
end

initval=A001_Images_Set_Experiment(exp);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname{ColourIdx}); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------

Nbac=Processing_Measure_Database(S)

[~,chan_no]=size(S) 

%[chan_no,~]=size(S)%%I THINK THIS WAS AN ERROR
thisisSizeS = size(S)
divtimes=[];
reptimes=[];
for ch=1:chan_no  %for each channel
 
chan_no-ch;
BF=S(ch).channels.kymo_BF;
BW=MakeBinaryEdgeImagefrom(BF, initval.BW_Threshold);
BW=BW';

if actions.plot %nargin<1,
    figure;
    subplot(1,2,1); P_Color(BF',500,1000, 'grey');;
    subplot(1,2,2); P_Color(1-BW,500,1000, 'grey');
   % [~]=ginput(1);
   pause(2);
%     demoloadpth='D:\jkerssemakers\My Documents\BN_All_ActiveProjects\BN_All_Programming\BN_All_CommonMatlab\JK_Matlab CleanTools';
%     demoloadnm='\DemoData_Bacteria';
%     save(strcat(demoloadpth,demoloadnm),'kymoBF_demo');    
end


[~,repno]=size(S(ch).channels.ReplicationCluster);
%clean the structure if needed
 if isfield(S(ch).channels, 'AutoDivision')
    S(ch).channels = rmfield(S(ch).channels, 'AutoDivision');
 end
 
for m=1:repno  %1) For every replicationcluster....
   %First, for convenience, store labels and existence of parent and
   %daughter cells
    ThisDivInfo=Add_FamilyMembers(S,M,ch,m);
    %----------------------------------------------------------------------
    %(if found) create 'left'and 'right' forks that start at termination and extend to the new
     %daugther init points. The new division should lie in the region encompassed by this fork and the subsequent region.
    %this division marks the end of this bacterium.   
    ThisDivInfo=Add_Forks(S,M,ch,ThisDivInfo);
    %Next, find Division time in an area aftr termination, bounded by the cluster positions of the daughter cels. 
    %The first detection of an edge count as division time. If nothing is
    %detected for the duration of the daughtercell's life, or no daughter
    %cells are found in the database, the division time is set as the
    %termination time
    ThisDivInfo=Add_division_time(S,M,ch,ThisDivInfo,BW);
    
    %Finally, add the division results to the database
   
    S(ch).channels.AutoDivision(m).family=ThisDivInfo.family;  
    S(ch).channels.AutoDivision(m).Forks=ThisDivInfo.Forks;  
    S(ch).channels.AutoDivision(m).divtype=ThisDivInfo.divtype; 
    S(ch).channels.AutoDivision(m).divisiontime=ThisDivInfo.divisiontime;
    S(ch).channels.AutoDivision(m).accepted=1; % just to be sure field exists
    S(ch).channels.ReplicationCluster(m).Forks=ThisDivInfo.Forks;

    %pcolor(1-BW); colormap grey, shading flat
    title(strcat('exp',exp,'channel' ,num2str(ch)));
    %[~]=ginput(1);
end

%_______________________________________________________________________
% Next, with all division times stored, we re-loop through the database 
%to find the birth times. 
%This is simply the division time of the parent.
%     If no data is stored of the parent, it is the initiation time
[~,repno]=size(S(ch).channels.ReplicationCluster);
for m=1:repno  %1) For every replicationcluster....
    ThisDiv=S(ch).channels.AutoDivision(m);
    ThisDivNw=Add_birth_time(S,M,ch,ThisDiv);    
    %Finally, add to database
    S(ch).channels.AutoDivision(m).birthtime=ThisDivNw.birthtime;  
    S(ch).channels.AutoDivision(m).birthtype=ThisDivNw.birthtype;
end
%------------------------------------------------------------------------

for m=1:repno  %1) For every replicationcluster....
    ThisDiv=S(ch).channels.AutoDivision(m);
    ThisDiv=Add_division_edges(S,M,ch,ThisDiv,BW);       
    %Finally, add to database
    S(ch).channels.AutoDivision(m).edges=ThisDiv.edges;
end


%plot & save section------------------------------------------------
%--------------------------------------------------------------------
%Show the cycle time as difference between birth and division
for m=1:repno  %1) For every replicationcluster.... 
   strt=S(ch).channels.AutoDivision(m).birthtype;
   stp=S(ch).channels.AutoDivision(m).divtype;
   ok=((strcmp(strt, 'OK') & strcmp(stp, 'OK')) &...
        S(ch).channels.AutoDivision(m).edges.edgesok);  %only the good ones!
   %ok=1;
   if ok
   d=S(ch).channels.AutoDivision(m).divisiontime-S(ch).channels.AutoDivision(m).birthtime;
   r=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(end)-S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(1);
   divtimes=[divtimes d];
   reptimes=[reptimes r];
   end
end



%Finally, save stuff
if actions.savedata
    save(outname,'S', '-append');
end
display('done');
end
[good,bad]=Processing_Measure_GoodDatabase(S)


figure;
subplot(1,2,1);
plot(reptimes,'o'); hold on;
plot(divtimes,'ro'); hold on;
%legend(['Replication', 'Automated Division Meas']);
xlabel('cell');
ylabel('time, frames');
subplot(1,2,2);
plot(reptimes,divtimes, 'k*'); hold on;
plot([0 110], [0 110], 'k-');
xlabel('replication time, frames');
ylabel('division time, frames');
end


function ThisDivInfo=Add_Forks(S,M,ch,ThisDivInfo)
%Create artifical 'fork' traces spannning the diffusion phase; to be used
%for setting domain borders for allocating division edges
    m=ThisDivInfo.family.me.idx;
    pp=ThisDivInfo.family.parent.idx;
    dL=ThisDivInfo.family.rightdaughter.idx;
    dR=ThisDivInfo.family.leftdaughter.idx;

    frs_m=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames;
    xs_m= S(ch).channels.ReplicationCluster(m).PosKyTracCom.trackpos;
    plot(frs_m,xs_m, 'o-'); hold on;
    fr_MT=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(end);   %mother termination frame
    x_MT=xs_m(end);  %mother termination pos
    if dL>0
       fr_LDI= S(ch).channels.ReplicationCluster(dL).PosKyTracCom.frames(1);     %Left daughter initiation frame 
       x_LDI=round(S(ch).channels.ReplicationCluster(dL).PosKyTracCom.trackpos(1));  % Left daughter initiation frame      
       frs_Lfork=fr_MT+1:fr_LDI-1;
       if length(frs_Lfork)==0,frs_Lfork=fr_MT;, end;      
       lf=length(frs_Lfork);     
       xs_Lfork=linspace(x_MT,x_LDI,lf);  
       Forks.left=[frs_Lfork' xs_Lfork' ];
       if lf>0
       plot(frs_Lfork,xs_Lfork,'ro-'); hold on;
       end
    else
        Forks.left=[fr_MT+1  x_MT];
        plot(fr_MT+1,  x_MT,'r*-'); hold on;
    end
   
    if dR>0
       fr_RDI= S(ch).channels.ReplicationCluster(dR).PosKyTracCom.frames(1);   %Right daughter initiation frame 
       x_RDI=round(S(ch).channels.ReplicationCluster(dR).PosKyTracCom.trackpos(1)); %Right daughter initiation frame
       frs_Rfork=fr_MT+1:fr_RDI-1;       
       if length(frs_Rfork)==0,frs_Rfork=fr_MT;, end;        
       lf=length(frs_Rfork);
       xs_Rfork=linspace(x_MT,x_RDI,lf);
       Forks.right=[frs_Rfork' xs_Rfork' ];  
       if lf>0
       plot(frs_Rfork,xs_Rfork,'mo-'); hold on;
       end
    else
        Forks.right=[fr_MT+1  x_MT];
        plot(fr_MT+1,  x_MT,'mo-'); hold on;
    end
    ThisDivInfo.Forks=Forks;
end   

function filim=Enhance_FilamentousStructures(KBF,tol);
    [r,c]=size(KBF);
    points=[];
    for i=1:r
        prf=KBF(i,:);
        mxi=F110_Get1DpeaksFlatBottom(prf,tol); %get peaks per line
        BW(i,mxi)=1;
        if length(mxi)>0
            ptsi=zeros(length(mxi),3);
            ptsi(:,1)=mxi;
            ptsi(:,2)=i;
            ptsi(:,3)=1;
            points=[points; ptsi];
        end
        %BW(i,2)=1;   %just to have points in every 'frame'
        %BW(i,end-2)=1;
        prf=BW(i,:);
        lm=length(mxi); 
    end
     fil_initval.GaussLongAxis=7;
     fil_initval.GaussShortAxis=3;
     fil_initval.GaussSearchRadius=8;
     fil_initval.weighpointsbyintensity=1;     
     filim=JKD2_IM_filament_it(points,KBF,fil_initval);
     filim=filim';
     if 0
      subplot(1,2,1); pcolor(KBF); shading flat; colormap hot
      subplot(1,2,2); pcolor(filim); shading flat; colormap hot
      [~]=ginput(1);
     end
end

function BW=MakeBinaryEdgeImagefrom(KBF,tol);
% This function transforms an image into an edge-detected, binary image;
% this is useful if peaks with strongly varying intensities exist. JK13

KBF=Enhance_FilamentousStructures(KBF,tol);
    
BW=0*KBF;
[r,c]=size(KBF);

points=[];
for i=1:r
    prf=KBF(i,:);
    mxi=F110_Get1DpeaksFlatBottom(prf,tol); %get peaks per line
    BW(i,mxi)=1;
    if length(mxi)>0
        ptsi=zeros(length(mxi),3);
        ptsi(:,1)=mxi;
        ptsi(:,2)=i;
        ptsi(:,3)=1;
        points=[points; ptsi];
    end
    %BW(i,2)=1;   %just to have points in every 'frame'
    %BW(i,end-2)=1;
    prf=BW(i,:);
    lm=length(mxi); 
end

  
   subplot(1,2,2); P_Color(1-BW,500,1000, 'grey');
  %Series of binary operations
  BW=bwmorph(BW,'clean');      
  BW=bwmorph(BW,'dilate',1);
                %BW=bwmorph(BW,'clean');
                %BW=bwmorph(BW,'erode',1);
  BW=bwmorph(BW,'skel', Inf); 

%     P_Color(1-BW,500,1000, 'grey');
% [~]=ginput(1);
end

    function ThisDivInfo=Add_division_time(S,M,ch,ThisDivInfo,BW)
    %Find Division time in an area aftr termination, bounded by the cluster positions of the daughter cels. 
    %The first detection of an edge count as division time. If nothing is
    %detected for the duration of the daughtercell's life, or no daughter
    %cells are found in the database, the division time is set as the
    %termination time

    m=ThisDivInfo.family.me.idx;
    pp=ThisDivInfo.family.parent.idx;
    dL=ThisDivInfo.family.rightdaughter.idx;
    dR=ThisDivInfo.family.leftdaughter.idx;
    
    ThisDivInfo.divtype='DNF';  %'did not finish'
    ThisDivInfo.divisiontime=0;
       
    replistart= S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(1);    %start time for searching
    replistop= S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(1);    %start time for searching
    
    Forks=ThisDivInfo.Forks;
    divtime=Forks.right(1,1); 
    if dR>0 && dL>0
    %f0= S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(1);    %start time for searching
    fL1= S(ch).channels.ReplicationCluster(dL).PosKyTracCom.frames(end); 
    fR1= S(ch).channels.ReplicationCluster(dR).PosKyTracCom.frames(end);

    %1 find proper frametime index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %make a left domain border (division should take place within these
    %borders)
    Lforkpos=Forks.left(:,2);
    Lforkfr= Forks.left(:,1);
    LtracePos=S(ch).channels.ReplicationCluster(dL).PosKyTracCom.trackpos;
    LtraceFr=S(ch).channels.ReplicationCluster(dL).PosKyTracCom.frames';
    leftborderpos=[Lforkpos ; LtracePos];
    leftborderfr=[Lforkfr ; LtraceFr];
    plot(leftborderfr,leftborderpos);
    %make a right domain border (division should take place within these
    %borders)
    Rforkpos=Forks.right(:,2);
    Rforkfr=Forks.right(:,1);
    RtracePos=S(ch).channels.ReplicationCluster(dR).PosKyTracCom.trackpos;
    RtraceFr=S(ch).channels.ReplicationCluster(dR).PosKyTracCom.frames';
    rightborderpos=[Rforkpos ; RtracePos];
    rightborderfr=[Rforkfr ; RtraceFr];
    plot(rightborderfr,rightborderpos);
    if 0
    plot(rightborderfr,rightborderpos, 'r-'); 
    pcolor(1-BW); colormap grey, shading flat; hold on
    [~]=ginput(1);
    end
    cnt=0; nodiv=1;
    divtime=Rforkfr(1);

    maxcnt=min(length(rightborderfr),length(leftborderfr));
    %start settings for division info
    divtime=Forks.right(1,1); 
    Rdivy=Forks.right(1,2);
    Ldivy=Forks.left(1,2);
    ThisDivInfo.divtype='DNF';  %'did not finish'
    SecCount=0;
    while nodiv  && cnt<maxcnt  %search until end of two daugthers   
    cnt=cnt+1;  
    [r,c]=size(BW);
    x=rightborderfr(cnt); 
    
    %Crop for minimum row indeximage;
    ylo=max(ceil(leftborderpos(cnt)), 1);
    yhi=max(floor(rightborderpos(cnt)),2);
    
    %Crop for maximum row indeximage;
     ylo=min(ylo, r-1);
     yhi=min(yhi,r);
     
    
    if yhi>ylo, trc=BW(ylo:yhi,x);  %select the binary image part encompassed by borders
    else trc=BW(ylo,x); end
%     if sum(trc)>0;
%         SecCount=SecCount+1;
%     end
    if sum(trc)>0 % && SecCount==10;
      nodiv=0;
      divtime=x; 
      Rdivy=ylo;
      Ldivy=yhi;
      ThisDivInfo.divtype='OK';    %'found a proper division'
%       plot([divtime divtime], [Ldivy Rdivy], 'r*-'); 
%       hold on;
     end
    end    
    end 
    ThisDivInfo.divisiontime=divtime;
    end    
    
    
    function ThisDivInfo=Add_FamilyMembers(S,M,ch,m);
  %%1 some genealogy: find names of the two daughters, and the parent ;
    %find if they exist
    [~,repno]=size(S(ch).channels.ReplicationCluster);
       
    ThisDivInfo.family.me.name=M(ch).channels.RepClicks(m).name;  
    ThisDivInfo.family.me.idx=m; 
    
    ThisDivInfo.family.parent.name=floor(ThisDivInfo.family.me.name/2);  %to find it in the database  
    ThisDivInfo.family.leftdaughter.name=ThisDivInfo.family.me.name*2+1;
    ThisDivInfo.family.rightdaughter.name=floor(ThisDivInfo.family.me.name*2);
    %---------------------------------------------------------------------

    ThisDivInfo.family.leftdaughter.idx=0;
    ThisDivInfo.family.rightdaughter.idx=0;
    ThisDivInfo.family.parent.idx=0;
    
    for j=1:repno 
        if M(ch).channels.RepClicks(j).name==ThisDivInfo.family.leftdaughter.name, 
            ThisDivInfo.family.leftdaughter.idx=j;
            dL=j; 
        end
        if M(ch).channels.RepClicks(j).name==ThisDivInfo.family.rightdaughter.name, 
            ThisDivInfo.family.rightdaughter.idx=j;
            dR=j; 
        end
        if M(ch).channels.RepClicks(j).name==ThisDivInfo.family.parent.name,
            ThisDivInfo.family.parent.idx=j;
            pp=j; 
        end
    end    
    end 
    
    
    function ThisDiv=Add_birth_time(S,M,ch,ThisDiv);
    % Next, we find the birth time. This is simply the division time of the parent.
    % If no data is stored of the parent, it is the initiation time    
    m=ThisDiv.family.me.idx;
    pp=ThisDiv.family.parent.idx;
    if pp>0 %if parent exists
         ThisDiv.birthtime= S(ch).channels.AutoDivision(pp).divisiontime;
         ThisDiv.birthtype='OK';
    else
        ThisDiv.birthtime=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(1);
        %This is just the initation time
        ThisDiv.birthtype='DNS';  %RDL changed from 'DNS' = 'did not start' to 'OK'
    end
    end
    
    function ThisDiv=Add_division_edges(S,M,ch,ThisDiv,BW);
    %This function  per cluster:    
    %a) find nearest edges left and right of the fluorescent traces consisting of
    %the cluster, its parent fork and the daughters, and between the start
    %time and the division time  
    %b) removes found edge points that are too far away, or can be otherwise
    %considered outliers
    %
    %cleaned edge points output is added as discontinuous traces to the divison data base. 
    %---------------------JacobKers 2013-------------------------
    
    %first, get some family info
    m=ThisDiv.family.me.idx;
    pp=ThisDiv.family.parent.idx;
    dL=ThisDiv.family.rightdaughter.idx;
    dR=ThisDiv.family.leftdaughter.idx;
    
    if 2*ThisDiv.family.parent.name==ThisDiv.family.me.name,
        leftdaughter =1;
    else
        leftdaughter =0;
    end
         
    %------------------------------------------------
    %get a continous time-pos trace of relevant cluster sections: 
    %fork parent, present cluster, fork daughter, daughter cluster
    %this total section encompasses the times that division edges can exist
    sectionL1_pos=[];  sectionL1_time=[];
    sectionL2_pos=[];  sectionL2_time=[];
    sectionL3_pos=[];  sectionL3_time=[];
    sectionL4_pos=[];  sectionL4_time=[];    
    if pp>0 
        if leftdaughter
         sectionL1_time=S(ch).channels.ReplicationCluster(pp).Forks.left(:,1);
         sectionL1_pos=S(ch).channels.ReplicationCluster(pp).Forks.left(:,2);  
        else
         sectionL1_time=S(ch).channels.ReplicationCluster(pp).Forks.right(:,1);
         sectionL1_pos=S(ch).channels.ReplicationCluster(pp).Forks.right(:,2);  
        end
    end
         sectionL2_time=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames';  
         sectionL2_pos=S(ch).channels.ReplicationCluster(m).PosKyTracCom.trackpos;
         sectionL3_time=S(ch).channels.ReplicationCluster(m).Forks.left(:,1);
         sectionL3_pos=S(ch).channels.ReplicationCluster(m).Forks.left(:,2);  
    
    if dL>0 
     sectionL4_time=S(ch).channels.ReplicationCluster(dL).PosKyTracCom.frames';
     sectionL4_pos=S(ch).channels.ReplicationCluster(dL).PosKyTracCom.trackpos;         
    end
   %-----------------------------------------------------------------------    
    %then, build a clusterline, consisting of
    %the cluster, its parent fork and the daughters, and between the start
    %time and the division time 
    Lclusters_pos=[sectionL1_pos' sectionL2_pos' sectionL3_pos' sectionL4_pos']';
    Lclusters_time=[sectionL1_time' sectionL2_time' sectionL3_time' sectionL4_time']';
    
    idx1=find(Lclusters_time==ThisDiv.birthtime);
    idx2=find(Lclusters_time==ThisDiv.divisiontime-1);
    
    Lclusters_time= Lclusters_time(idx1:idx2);
    Lclusters_pos= Lclusters_pos(idx1:idx2);
   %-----------------------------------------------------------------------  
   
   %for every time point, search the nearest 'left' or 'right' point in the
   %binary edge image
   [kymlength,kymdur]=size(BW);
   Le=length(Lclusters_time);
   left=zeros(Le,1);
   right=zeros(Le,1);
   close all;
   for i=1:Le
       fr=min(Lclusters_time(i), kymdur);
       clusterpos=min(ceil(Lclusters_pos(i)),kymlength);
       clusterpos=max((clusterpos),1);
       %Get nearest left point
       BWleftline=BW(clusterpos:end,fr);
       sel=find(BWleftline==1);  %get peak indices
       if ~isempty(sel), left(i)=clusterpos+sel(1);   %pos of first nonzero element
       else              left(i)=kymlength; end
       %Get nearest right point 
       
       BWrightline=BW(1:clusterpos,fr);
       sel=find(BWrightline==1);  %get peak indices
       if ~isempty(sel), right(i)=sel(end);   %pos of first nonzero element
       else              right(i)=1; end
   end 
    
   %--------------------------------------------------------
    if 0 %plot black and white picture with one curve overlay
    figure;
    plot(Lclusters_time, Lclusters_pos, 'r-'); hold on; %search center line
    plot(sectionL2_time, sectionL2_pos, 'ro-'); hold on; %replication cluster
    plot(Lclusters_time, left, 'o-'); hold on; %left edge
    plot(Lclusters_time, right, 'mo-'); hold on; %right edge
    pcolor(1-BW); colormap grey, shading flat; hold on
    %[~]=ginput(1);
    end
    
    %--------------------------------------------------------------
    %In the following: remove found edge points that are too far away, or can be otherwise
    %considered outliers   
    
    %first, collect edges
    edges.repfrs=sectionL2_time;
    edges.reppos=sectionL2_pos;
    edges.frs=Lclusters_time;
    
    edges.left=left; 
    edges.right=right;
        
    %evaluate edge points on criteria such as minimal length of the
    %bacterium; expand the 'edge' properties with fits, acceptance flags
    %etc.   
    edges=Expand_Edges(edges,BW); 
       
    %finally, allocate:
    ThisDiv.edges=edges;
    end
    
        
    function edges=Expand_Edges(edges,BW);
    close all;  
    
    %----------------------------------------------------------------------
    % RdL -- Making bacterial lifetime longer than 
    % cycle for lengthened BacPics -- 
    
    Nfrs=0; % Number of frames added to cycle
    Adum=edges.frs; %create dummy vars
    Ldum=length(Adum);
    
    %redefine edges.frs with prolonged length
    edges.frs2=linspace(Adum(1),Adum(Ldum)+Nfrs,Adum(Ldum)+Nfrs-Adum(1)+1)';
    %----------------------------------------------------------------------
    
    lifelength=length(edges.frs2);
    edges.acceptpoints=zeros(lifelength,1);   %flags used to indicate if points are ok
    edges.edgesok=1;                   %flags used to indicate if bac edges are ok
    
    midbac=(edges.left+edges.right)/2;   %midline of bacterium
    baclength=edges.left-edges.right;    %length of bacterium     
    
    %First, find and reject too short lengthts per time point
    sel=find(baclength>7); edges.acceptpoints(sel)=1;
     
    %Next, fits on the remaining 'good' points
    Lpars=polyfit(edges.frs2(sel),edges.left(sel)-midbac(sel),1); %linear fit on growth
    Rpars=polyfit(edges.frs2(sel),edges.right(sel)-midbac(sel),1); %linear fit on growth (just the negative)
    Mpars=polyfit(edges.frs2(sel),midbac(sel),2);  %polynomial fit on absolute position
    
    edges.midfit=Mpars(1)*edges.frs.^2+Mpars(2)*edges.frs+Mpars(3);
    edges.leftfit=Lpars(1)*edges.frs+Lpars(2)+edges.midfit;
    edges.rightfit=Rpars(1)*edges.frs+Rpars(2)+edges.midfit;
    
    edges.midfit2=Mpars(1)*edges.frs2.^2+Mpars(2)*edges.frs2+Mpars(3);
    edges.leftfit2=Lpars(1)*edges.frs2+Lpars(2)+edges.midfit2;
    edges.rightfit2=Rpars(1)*edges.frs2+Rpars(2)+edges.midfit2;
    
    %Parameters of edges
    edges.Lpars=Lpars;
    edges.Rpars=Rpars;
    edges.Mpars=Mpars;
    
    %Based on the fits on the 'good' points, evaluate this bacterial cycle
    if Lpars(1)<0, 
        edges.edgesok=0; end  %We reject shrinking bacteria
    if max(edges.midfit)-min(edges.midfit)>75,  
        edges.edgesok=0; end  %We reject bacteria wandering too much 
    if min(edges.leftfit-edges.rightfit)<2, 
        edges.edgesok=0; end  %We reject weird growth beahivior
    if lifelength>100,
        edges.edgesok=0; end  %We endless existence
    
    if 0 %show the results; relative pane
    plot(edges.frs,edges.left-midbac,'ro-'); hold on; %all points
    plot(edges.frs,edges.right-midbac,'ro-');
    plot(edges.frs(sel),edges.left(sel)-midbac(sel),'o-'); hold on; %good points
    plot(edges.frs(sel),edges.right(sel)-midbac(sel),'o-');
    plot(edges.frs,edges.leftfit-edges.midfit,'k-'); hold on;             %fit
    plot(edges.frs,edges.rightfit-edges.midfit,'k-');
    [~]=ginput(1);
    end
    
    if 0 %show the results; absolute pane
    hold on
    plot(edges.frs,edges.leftfit,'r-');  %fit
    plot(edges.frs,edges.rightfit,'r-'); 
    pcolor(1-BW); colormap grey, shading flat; hold off
    [~]=ginput(1);   
    end
    
    if 0 %plot relative -- RdL
        hold on
        plot(edges.midfit2)
        plot(edges.leftfit2)
        plot(edges.rightfit2)
        plot(edges.left)
        plot(edges.right)
        plot(edges.reppos)
        hold off
    end
    edges=edges;
    end
