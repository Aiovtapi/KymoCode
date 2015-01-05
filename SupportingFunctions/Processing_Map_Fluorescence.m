function dum=Processing_Map_Fluorescence(ThisBac,ThisRep,stripmov_FL,kymo_FL,initval);
%Function makes a map of fluorescence and division props

%frames and positions------------------------------------------------------
frs=ThisBac.PosKyTrac.frames;      
lft=ThisBac.PosKyTrac.left;
rht=ThisBac.PosKyTrac.right;
repfrs=ThisRep.PosKyTracCom.frames;

%some more props
middiv=(rht+lft)/2;  %midline bactrium

%select a time slot encompassing replication and division cycle
%--------------------------------------------------------------------------
xt=initval.extension;  %extension before and after replication and division times
hifr=repfrs(end)+xt;
lofr=repfrs(1)-xt;
frs_ext=[lofr:hifr];

%fits to midline, second order---------------------------------------------
ppM=polyfit(frs,middiv,2);
fit_mid=ppM(1)*(frs_ext).^2+ppM(2)*frs_ext+ppM(3); 

%fits RELATIVE to MIDLINE-----------------------
ppL=polyfit(frs,lft-middiv,1);
%         fit_lft=ppL(1)*(frs_ext).^2+ppL(2)*frs_ext+ppL(3);    
fit_lft=round(ppL(1)*(frs_ext)+ppL(2));  

ppR=polyfit(frs,rht-middiv,1);
%fit_rht=ppR(1)*(frs_ext).^2+ppR(2)*frs_ext+ppR(3);       
fit_rht=round(ppR(1)*(frs_ext)+ppR(2)); 

%using fit results, we now fetch a fluorescence profile from the kymograph:
%----------------------------------------------
hor=max(fit_rht-fit_lft)+4;
ver=hifr-lofr+1;
pic=zeros(ver,hor);

[kr,kc]=size(kymo_FL);
for t=1:ver
fri=max(frs_ext(t),1); fri=min(fri,kr);  %inbound frame number
loix=max(fit_lft(t)+fit_mid(t), 1);  %inbound low pos index
hiix=min(fit_rht(t)+fit_mid(t), kc); %inbound hi pos index
idx=round([loix:hiix]);% indices for kymograph line
prof=kymo_FL(fri,idx); prof=prof-min(prof)+0.001;
lixx=length(idx); shft=(hor-lixx)/2;
idx2=round([1:lixx]+shft); %centered indices;
pic(t,idx2)=prof; 

%----------------------------------------------------------
end

%plot menu;  ---------------------------

pcolor(pic'); colormap hot; shading flat;  hold on
plot(frs-frs_ext(1)+1,rht-middiv+hor/2+1.5,'-o'); 
plot(frs-frs_ext(1)+1,lft-middiv+hor/2+1.5,'-o');

title(strcat('Bacterium'));
xlabel('time, frames');
ylabel('position, pixels');
axis([1 ver 1 hor]);
%[~]=ginput(1);
dum=1;
