%function dum=Division00_Analyze_Kymo(kymo_BF,initval)
 kymo=kymo_BF2 ; kymo=kymo(1:550,:);
 kymo_FL=kymo_FL; kymo_FL=kymo_FL(1:550,:);

close all; 
figure; 
subplot(1,3,1);pcolor(kymo_FL); shading flat; colormap hot; title('Fluorescence Kymograph');
subplot(1,3,2); pcolor(kymo); shading flat; colormap hot; title('BrightField Kymograph');
 % figure; plot(kymo(1,:),'-o')

bacprops.count=8;    %initial number in channel (rounded to upper 2^n (used for labeling)
bacprops.exitpos=c-2;
bacprops.length=20; 
bacprops.sim=0;
bacprops.bacgro=0.15;
bacprops.whenexit='left' %'right'
[r,c]=size(kymo); 
 
 
Division=Division_Init(bacprops);

bacprops.treshold=Get_edge_treshold_kymo(kymo);
bacprops.FLtreshold=Get_edge_treshold_kymo(kymo_FL);

bacprops
for fr=2:r;
    fr;
    bacprops.count;
%per frame:
    %1)find existing bacteria and store some current properties (like
    %number in channel etc.)
        framestate=Division_GetLivingEcoli(Division);   %Get properties of current bacteria in channel (from current database)  
        framestate.frame=fr;
    %per bacterium:
        for i=1: framestate.livingbacno
            bacno=framestate.livingbac(i);   %pick one of the 'living' bacteria        
            Division=Division_UpdateBacPos(Division,framestate,bacprops,bacno,kymo);    %update its position
            Division=Division_CheckExit(Division,framestate,bacprops,bacno); %Check if exit was reached
            Division=Division_CheckEnd(Division,framestate,bacno,kymo); %Check if end of measurement was reached          
           [Division,bacprops]=Division_Checkdiv(Division,framestate,bacprops,bacno,kymo); %Check if division took place 
           
           %(this will end its 'live' and initiate its offspring)        
           %[Division,bacprops]=CheckFluo(Division,framestate,bacprops,bacno,kymo_FL); %Check when spots are present
           
        end
end
Bax = struct2cell(Division);
fate=squeeze(Bax(5,:,:));   
subset=find(strcmp(fate, 'nonexistent')~=1);  %keep indices of present bacteria
Division=Division(subset);
[le,~]=size(Division);
figure;



%Further analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:le
    %bacno=subset(i);
   % subplot(1,3,3);  
    bacno=i;
    life=length(Division(bacno).left); 
    ax=Division(bacno).frames;
    sticks=[Division(bacno).left Division(bacno).right]';
    switch mod(Division(bacno).name,2)
        case 0, clr='r-*';
        case 1,clr='g-*';
    end
    ax2=[ax ax]';
    plot(sticks,ax2, clr); hold on
    plot(sticks,ax2, 'k*'); hold on
    axis([0 c 0 fr]);
    title('Edges Analysis');
    
end




% %check replication time
% Bax = struct2cell(Division);
% fate=squeeze(Bax(5,:,:));   
% subset=find(strcmp(fate, 'divided')==1);  %find indices of present bacteria
% Division=Division(subset);
% 
% ls=length(subset);
% spotlife=[];
% divtime=[];
% divlength=[];
% for i=1:ls
%     bacno=i;
%     t=length(Division(bacno).left);
%     divtime(i)=t;
%     gro=Division(bacno).right-Division(bacno).left;
%     divlength(i)=gro(t);
%     
% end
% spotlife=spotlife*2.5; %minutes
% divlength=divlength*0.16; %minutes
% divtime=divtime*2.5 % minutes
% 
% figure;
% subplot(2,2,1);
% plot(divtime, divlength, 'o');
% title('end size vs. division time');
% xlabel('division time,minutes'); ylabel('end size, microns');
% 
% subplot(2,2,3);
% binz=6; 
% %hx=(min(divtime):(range(divtime))/binz:max(divtime));   %make an axis
% hx=(0:200/binz:200);   %make an axis
%     sthst=hist(divtime,hx);
%     bar(hx,sthst);
%     title('division time'); 
%     xlabel('division time, minutes'); ylabel('counts');
% 
%     
%     
%     
% 
% subplot(2,2,4);
% binz=6;
% hx=(0:200/binz:200);   %make an axis
% 
%     sthst=hist(spotlife,hx);
%     bar(hx,sthst);
%     title('replication time'); 
%     xlabel('replication time, minutes'); ylabel('counts');
% 
% 
% 
% %for info: Some basic database actions %%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Bax = struct2cell(Division);
% TokenBac=squeeze(Bax(:,2,:));   %Example how to pick one bacterium state
% fate=squeeze(Bax(6,:,:));   %example how to pick one property
% subset=find(strcmp(fate, 'alive'));  %find indices of present bacteria
% Division(subset).name ;                    %...and their names
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%