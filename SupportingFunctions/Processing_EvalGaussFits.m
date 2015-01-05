function [f1out]=Processing_EvalGaussFits(f1,FLspot, fluofull,initval);
%evaluating 1D Gauss peaks 
%We judge the peak pairs based on what we know of the replication cycle
% -spots will not be outside cell
% -spots will be strong enough
% -spots will be equal enough

%Input:
%FLspot: image with spots 


% fluofull = 
% 
%     area_bac: 64
%     area_spot: 57
%     content_cytoplasm1: 48897.5239552546
%     content_cytoplasm2: 78485.1472936904
%     content_signal: 277678.914334679
%     content_spots1: 228781.390379424
%     content_spots2: 199193.767040988
%     content_total: 386145.718331106
%     curve_medianofmax: [1x23 double]
%     curve_medianofsum: [1x23 double]
%     level_dark: 154.701802680748
%     level_fluotreshold: 402.450588337163
%     level_medianofmax: 915.226862192169
%     level_medianofmax_yposcurve: [1x23 double]
%     level_medianofsum: 3363.09749257291
%     level_peak: 17130.1088949146
%     noise_dark: 123.874392828207
%     peak_xpos: [1x23 double]
%     peak_ypos: [1x23 double]
%     ratio_FS: 0.71735287325744
%     ratio_SN: 13.6947337582286
%     wherebac: [64x1 double]
%     wheredark: [286x1 double]
%     wherefluo: [151x1 double]
%     wherespot: [57x1 double]


%output:
% f1 out = [X0,X1,Y0,Y1, Background amplitude,Peak0, Peak1 spotflag1 spotflag 2] 
%with 'spotflag' =1 if this is a good spot, otherwisew 'NaN')

%f_out:= xmainpeak xsecondarypeak ymainpeak ysecondarypeak background
f1out=f1;
f1out.Spot1OK=1;
f1out.Spot2OK=1;
[r,c]=size(FLspot);

%1) similar amplitudes? ----------------
if f1.contentspot1/f1.contentspot2<0.2, f1out.Spot1OK= NaN;end
if f1.contentspot2/f1.contentspot1<0.2, f1out.Spot2OK= NaN;end


%spots in range?
if ~(f1.X0>1)*(f1.X0<c)*(f1.Y0>1) *(f1.Y0<r), f1out.Spot1OK= NaN;end
if ~(f1.X1>1)*(f1.X1<c)*(f1.Y1>1) *(f1.Y1<r), f1out.Spot2OK= NaN;end

%spots not too low -or negative- anyway?
if f1.contentspot1/fluofull.content_signal<0.05, f1out.Spot1OK= NaN;end
if f1.contentspot2/fluofull.content_signal<0.05, f1out.Spot2OK= NaN;end
dum=1;

