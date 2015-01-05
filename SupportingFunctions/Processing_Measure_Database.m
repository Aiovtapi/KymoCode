function Nbac=Processing_Measure_Database(S)
%measure size of complete database, given completeness of sets
Nbac=0;
[~,chan_no]=size(S);    
for i=1:chan_no  %for each channel
Rep=S(i).channels.ReplicationCluster;
[~,repno]=size(Rep);
for j=1:repno  %for each bacterium
Nbac=Nbac+1; end %set includes linked replication cycle
end
end