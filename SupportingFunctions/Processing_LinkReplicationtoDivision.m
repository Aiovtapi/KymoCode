function S=Processing_LinkReplicationtoDivision(S);
%Find LinkedReplicationCycle

%Link division to replication cycles
%rule: there is one bacterium for which a particular replication cycle is
%alone for a certain amount of time and vice versa. We'd like to link these
%two and add this to a new field in the 'Division' Data base


showplot=0;

[~,chan_no]=size(S);
for i=1:chan_no  %for each channel
chan_no-i
Div=S(i).channels.Division;                 %division data per channel
Rep=S(i).channels.ReplicationCluster;       %replication data per channel

[~,bacno]=size(Div);
[~,repno]=size(Rep);
for j=1:bacno  %for each bacterium; load (position edges,time data)
%'bacno', bacno-j
spotno=[];
ThisBac=Div(j);

frs=ThisBac.PosKyTrac.frames; divtime=frs(end); birthtime=frs(1);     
lft=ThisBac.PosKyTrac.left;
rht=ThisBac.PosKyTrac.right;
%check for each replication cycle if it fits in this bacerium; 
%first collect the (possibly multiple)cycles:-------------
for k=1:repno     
ini=Rep(k).PosKyTracCom.frames(1);    %inititiation
ter=Rep(k).PosKyTracCom.frames(end);  %termination
%
cond1=(ini<divtime)&(ter>birthtime);  %first: is there any temporal overlap?
if cond1
reppos=Rep(k).PosKyTracCom.trackpos;
mrep=mean(reppos);
mL=mean(lft);
mR=mean(rht);
cond2=((mrep>mL)&(mrep<mR));          %second: is there any spatial overalp?
if cond2
spotno=[spotno; [k ini]];
end
end
end 
if ~isempty(spotno)
%cleanup
ls=length(spotno(:,1));
spotno2=[];
for m=1:ls
ix=spotno(m,1);
repfrs=Rep(ix).PosKyTracCom.frames;
reppos=Rep(ix).PosKyTracCom.trackpos;

inrange=Processing_Check_Inrange(repfrs,reppos,frs,lft,rht); %final point....
...by point check: replication spots must lie within this bacteriums edges

if inrange,spotno2=[spotno2 ; spotno(m,:)];end
end
if ~isempty(spotno2)
 %spot life with earliest occurence is linked to this bacterium 
...(since we assume last termination is always equal or before first division:
[val,idx]=min(spotno2(:,2)); 
fitrep=spotno2(idx,1);
else
fitrep=0;  %no match found
end
else fitrep=0;  %no spot found
end
S(i).channels.Division(j).linkedrep=fitrep;
if fitrep>0
S(i).channels.ReplicationCluster(fitrep).linkedbac=j;
Rep(fitrep).linkedbac=j;
end

Div(j).linkedrep=fitrep;

%Extra action: by default at this point, all bacteria are accepted
S(i).channels.Division(j).accepted=1;
S(i).channels.DivClicks(j).accepted=1;


%-------optional plot menu
if showplot
if fitrep>0 
close(gcf);
FittedRep=Rep(fitrep);
plot(frs,rht,'k-*'); hold on
plot(frs,lft,'k-*');
repfrs=FittedRep. PosKyTracCom.frames;
reppos=FittedRep. PosKyTracCom.trackpos;
plot(repfrs,reppos,'-r*');
[~]=ginput(1);
end
%-------------------
end
end    
end


