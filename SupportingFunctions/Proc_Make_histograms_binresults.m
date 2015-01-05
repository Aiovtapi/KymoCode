function dum=Proc_Make_histograms_binresults(binresults,slots);
%Make histograms of binned results. JacobKers 2013
%---------------------------------------------------
%In this section, we make  stacked histograms of the data points over a period of
%normalized time, for example from t=0.2 to 0.8, where t=0 indicates
%initiation and t=1 indicates termination. The aim is to see if outliers
%affect our average & sigma too much. Note : this action only makes sense
%if we use relative times 
%example slots=[0.1  0.5  0.9];  %these are begin and end points of the relative time intervals
%Contents of 'results'
% binresults.losig;          %one sigma below the average, per bin
% binresults.av;             % the average, per bin
% binresults.hisig;          %one sigma above the average, per bin
% binresults.binaxis;   %values of the non-empty bins, per bin
% binresults.binaxis_all; %value of  bin, per point
% binresults.binaxis_all_scat; %random number to spread bin points (for plotting purposes)
% binresults.binvalues_all; %values themselves
% binresults.bincounts=traceno;

if nargin<2  %TEST modus with 50 fake cels
    close all; 
    slots=[0.1 0.3 0.6 0.9];
    binresults.binaxis_all=repmat(linspace(-0.1,  1.1,100),1,50);
    binresults.binvalues_all=repmat(3+randn(1,100),1,50);
end

figure;
binz=40;
lobin=0;
%hibin=max(binresults.binvalues_all);
hibin=200;
valaxis=linspace(lobin,hibin,binz);
slotno=length(slots)-1;
for i=1:slotno
    sel=find(binresults.binaxis_all>slots(i) & binresults.binaxis_all<slots(i+1));
    buf=binresults.binvalues_all(sel);
    hist_slot=hist(buf,valaxis);
    subplot(slotno,1,i);
    bar(valaxis,hist_slot); hold on;
    text(lobin+0.6*(hibin-lobin),0.8*max(hist_slot),strcat('From',num2str(slots(i)),'to',num2str(slots(i+1))));
    if i==1, title('distribution per relative time {init-ter}'), end
    if i==slotno; xlabel('Label no., a.u.'),    end
    axis tight
end
dum=1;

