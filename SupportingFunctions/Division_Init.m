function [Division,DivClicks]=Division_Init(kymo_BF,kymoprops,initval);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allocate first bacteria (to be replaced by detection, later)
% •	S(1).channels.Division   entries
% o	name: genealogy label; mother cells are2^n 
% o	fate: ‘exit’, ‘divided’ etc. used for removing incomplete cycles
% o	Linkedrep: database index of associated replication cycle (to be added
% later)
% o	PosClick: cycle start&finish times and positions, acquired by manual clicks
% o	PosKyTrac: tracked positions by analysing kymographs

Division=[];
DivClicks=[];

%later to be used in init)
n_tracks=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%single ['1010111'] style numer indicating family tree position
%number length-1= # of prior divisions
            %nth digit: left(1)- or right(0) division in nth division
            %(channel exit is on the far right)

%set names, start positions (yes, rather explicitly)
left=[];
right=[];
div=[];
stop=0;
titl=strcat('Initialization: Click Start positions, then right-click channel end position');
P_Color(kymo_BF(1:initval.startzoom,:),initval.kymolength,initval.startzoom,'jet'); title(titl); hold on;

while stop==0
    [x_L,y_L,but]=ginput(1);
    plot(x_L,y_L,'wo','MarkerSize', 12,'MarkerFace', 'r'); hold on;
    left=[left x_L];
    [x_R,y_R,but]=ginput(1);
    plot(x_R,y_R,'ro','MarkerSize', 10, 'MarkerFace', 'b'); hold on;
    right=[right x_R];
    plot([x_L x_R],[y_L y_R],'w-'); hold on;
    div=[div (y_L+y_R)/2];
    if x_L>kymoprops.width, stop=1; end
end
le1=length(left)-1;    %starting number, 
le2=2^(ceil(log2(le1))); %2^n to keep labelling consistent

padd_L=left(le1+1)*[1:le2-le1];
padd_R=right(le1+1)*[1:le2-le1];
paddif=1*[1:le2-le1];

left=ceil([left(1:le1) padd_L]);
right=ceil([right(1:le1) padd_R]);
div=ceil([div(1:le1) paddif]);

for i=1:le2-1

DivClicks(i).name=le2+i-1;
DivClicks(i).fate='alive';
DivClicks(i).linkedrep=0;

DivClicks(i).PosClick.firstfr=div(i);
DivClicks(i).PosClick.lastfr=div(i)+1;
DivClicks(i).PosClick.firstleft=left(i); 
DivClicks(i).PosClick.firstright=right(i);
DivClicks(i).PosClick.lastleft=left(i); 
DivClicks(i).PosClick.lastright=right(i);

Division(i).PosKyTrac.frames=div(i);
Division(i).PosKyTrac.left=left(i); 
Division(i).PosKyTrac.right=right(i);

end
for j=le1:le2, Division(j).fate='exit';  end     %(never seen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(gcf);
dum=1;
