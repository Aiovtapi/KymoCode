function twopoints=Adjust_twopoints(rowcor,colcor,twopoints,radz);
%apply correction values after proper rotation, rather spelled out...

%original positions
r=twopoints(1,1); 
c=twopoints(1,2);
r2=twopoints(2,1); 
c2=twopoints(2,2);

xcor=colcor*cos(radz)-rowcor*sin(radz);  %image horizontal correction 
ycor=colcor*sin(radz)+rowcor*cos(radz);  %image vertical correction 

c=c+xcor; 
c2=c2+xcor;
r=r+ycor; 
r2=r2+ycor;

twopoints(1,1)=r; 
twopoints(1,2)=c;
twopoints(2,1)=r2; 
twopoints(2,2)=c2;

dum=1;

