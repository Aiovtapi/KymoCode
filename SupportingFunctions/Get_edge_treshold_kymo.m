function treshold=Get_edge_treshold_kymo(kymo);
%This function detrmines a no-edge noise level for a kymograph
 [flag,cleandata]=Outlier_Flag(kymo,2,0.9,'all',0);
    treshold=3.5*std(cleandata);
end