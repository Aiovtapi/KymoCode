function [good, bad]=Processing_Measure_GoodDatabase(S)
%measure size of complete database, given completeness of sets
good.birth=0;
good.birthdiv=0; 
good.birthdivedges=0;
good.birthdivedgesuser=0;
bad.birth=0;
bad.birthdiv=0; 
bad.birthdivedges=0;
bad.birthdivedgesuser=0;
[~,chan_no]=size(S);    
for ch=1:chan_no  %for each channel
Div=S(ch).channels.AutoDivision;
[~,bacno]=size(Div);
for m=1:bacno  %for each bacterium   
    if ~isfield(S(ch).channels.AutoDivision(m), 'accepted')
        S(ch).channels.AutoDivision(m).accepted=1;  %To ensure field exists
    end
    
    ok1=strcmp(S(ch).channels.AutoDivision(m).birthtype, 'OK');    %birth ok
    ok2= strcmp(S(ch).channels.AutoDivision(m).divtype, 'OK');    ;%division ok
    ok3=S(ch).channels.AutoDivision(m).edges.edgesok;  %edges ok 
    ok4=S(ch).channels.AutoDivision(m).accepted;  %user accept


    if ok1,     good.birth=good.birth+1;  %set includes linked replication cycle
    else       bad.birth=bad.birth+1; end
    if ok1&ok2,     good.birthdiv=good.birthdiv+1;  %set includes linked replication cycle
    else       bad.birthdiv=bad.birthdiv+1; end
    if ok1&ok2&ok3,     good.birthdivedges=good.birthdivedges+1;  %set includes linked replication cycle
    else       bad.birthdivedges=bad.birthdivedges+1; end
    if ok1&ok2&ok3&ok4,     good.birthdivedgesuser=good.birthdivedgesuser+1;  %set includes linked replication cycle
    else       bad.birthdivedgesuser=bad.birthdivedgesuser+1; end
end
end