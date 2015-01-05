function dum=Processing00_KymoGraphTracker
%Automated Kymograph Analysis - Try U-track on this
%Jacob Kerssemakers,
%TNW-BN-ND lab 2012; developed for Charl Moolman


exp='TEST';

actions.loaddatabase=1; %default=1 (new analysis)

initval=A001_Images_Set_Experiment(exp);

%load the databases--------------------------------------------------
outname=strcat(initval.basepath,initval.outname); %processed inputs
outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
if actions.loaddatabase
load(outname,'S');
load(outname_usr,'M');
end
%------------------------------------------------------------------

%work through all 'accepted' bacteria cycles------------------------
[~,chan_no]=size(S);

pc=[];

for i=2:chan_no  %for each channel
KBF=S(i).channels.kymo_BF;
BW=0*KBF;
P_Color(KBF,500,500,'hot');
[~]=ginput(1);
[r,c]=size(KBF);
for i=1:r
    prf=KBF(i,:);
    mxi=F110_Get1DpeaksFlatBottom(prf,0.5); %get peaks per line
    BW(i,mxi)=1;
    BW(i,2)=1;   %just to have points in every 'frame'
    BW(i,end-2)=1;
    prf=BW(i,:);
    lm=length(mxi); 
end

  BW=bwmorph(BW,'clean');
  %movieInfo(i).xCoord = [peaks' 0*peaks'] ;
        
  BW=bwmorph(BW,'dilate',1);
  %BW=bwmorph(BW,'clean');
  %BW=bwmorph(BW,'erode',1);
  BW=bwmorph(BW,'skel', Inf);  
P_Color(1-BW,500,1000, 'grey');
[~]=ginput(1);


props.n_frames=r;       
initval.processeddatapath=cd;
%Make a  inputdataset compatible w/ utrack stuff------------------ 
  movieInfo = repmat(struct(...      
                'xCoord',[],....
                'yCoord',[],...
                'amp',[]),props.n_frames,1);
    for i=1:props.n_frames
    prf=BW(i,:);
    pks=find(prf==1);
        movieInfo(i).xCoord = [pks' 0*pks'] ;
        movieInfo(i).yCoord = [0*pks'  0*pks'];
        movieInfo(i).amp =    [0*pks'+1  0*pks'] ;    
    end
 %-----------------------------------------------------------------------
CustomTrack=Utrack_Shell(movieInfo,props,initval);
figure;
plot_Utrack_trajectories(CustomTrack);
[~]=ginput(1);
end
end




