function [spotno,f0out,f1out]=Processing_Judge_fitnumbers(f0,f1,im)
%Simple evaluating sorting of peaks 
%f_out:= xmainpeak xsecondarypeak ymainpeak ysecondarypeak background
[r,c]=size(im);

ampli_ratio=max([f1(6) f1(7)])/min([f1(6) f1(7)]);
diffu_ratio=max([f1(6) f1(7)])/(f1(5));

%similar amplitudes?  DISABLED
if ampli_ratio<5 
    cond.ampli=1 ;
else
    cond.ampli=0; 
end

%spot 1 in range?
if f1(1)>1 & f1(1)<c, 
    cond.inrange1=1; 
else
    cond.inrange1=0;
end

%spot 2 in range?
if f1(2)>1 & f1(2)<c
    cond.inrange2=1;
else
    cond.inrange2=0; 
end

%set spot number
if cond.inrange1*cond.inrange2*cond.ampli
    spotno=2; 
else
    spotno=1;
end

%background low enough? (to avoid diffuse sections w/ no prominent spots)
%DISABLED
%if diffu_ratio<3; spotno=0;end

%sort peaks
switch spotno
    case 2
    if f1(6)>f1(7)  %first peak brightest
    f0out=[f0(1) f0(2) f0(3) f0(4) f0(5) f0(6) f0(7)];
    f1out=[f1(1) f1(2) f1(3) f1(4) f1(5) f1(6) f1(7)];
    else  %swap peaks to make first peak brightest
    f0out=[f0(2) f0(1) f0(4) f0(3) f0(5) f0(7) f0(6)];
    f1out=[f1(2) f1(1) f1(4) f1(3) f1(5) f1(7) f1(6)];
    end
    case 1
    if cond.inrange1  %first peak in-range
    f0out=[f0(1) 0 f0(3) 0 f0(5) f0(6) 0];
    f1out=[f1(1) 0 f1(3) 0 f1(5) f1(6) 0];
    end  
    if cond.inrange2  %first peak largest&in-range
    f0out=[f0(2) 0 f0(4) 0 f0(5) f0(7) 0];
    f1out=[f1(2) 0 f1(4) 0 f1(5) f1(7) 0];
    end
    if ~cond.inrange1&~cond.inrange2
        f0out=[0 0 0 0 0 0 0 ];
        f1out=[0 0 0 0 0 0 0 ]; 
    end
    case 0
    f0out=[0 0 0 0 0 0 0 ];
    f1out=[0 0 0 0 0 0 0 ];
 end
    



