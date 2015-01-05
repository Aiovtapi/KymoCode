function ThisRepout=Processing_Get_InitTerm_Extensions(ThisRep,initval);
%This function analyzes the first part of a replication cycle in closer
%detail. JacobKers 2013

%fit settings user%--------------
hwix=9;   %spot-roi limits, vertical
hwiy=9;  %spot-roi limits, lateral
%-------------------------


xt=initval.extension;  %extension before and after replication and division times

frs_clc=ThisRep.PosKyTracCom.frames';  %frames between clicked positions
pos_clc=ThisRep.PosKyTracCom.trackpos;  %1D COM position s
pp=polyfit(frs_clc,pos_clc,1);

%Define time points before and after start and stop clicked positions
frs_ext_init=[frs_clc(1)-xt:frs_clc(1)-1]'; 
frs_ext_term=[frs_clc(end)+1:frs_clc(end)+xt]'; 

%back-extrapolation of first position, based on average trend--------------
posfit_ext_init=pp(1)*(frs_ext_init)-pp(1)*(frs_ext_init(end))+pos_clc(1); 

%forward-extrapolation of last position, based on average trend-----------
posfit_ext_term=pp(1)*(frs_ext_term)-pp(1)*(frs_ext_term(end))+pos_clc(end); 

%combine three sections in one encompassing trace--------------------------
frs_ext=[frs_ext_init' frs_clc' frs_ext_term']';
pos_ext=[posfit_ext_init' pos_clc' posfit_ext_term']';

%make sure negative times are avoided--------------------------------------
sel=find(frs_ext<1);
buf=find(frs_ext==1);
frs_ext(sel)=frs_ext(buf); 
pos_ext(sel)=pos_ext(buf);

%plot(frs_ext,pos_ext, '-o');  hold;

%store extended data separately
ThisRepout.PosKyTracCom.frames_ext=frs_ext;  %frames between clicked positions
ThisRepout.PosKyTracCom.trackpos_ext=pos_ext;  %1D COM positions

dum=1;