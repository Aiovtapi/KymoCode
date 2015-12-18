function [aa,ff,drift]=Get_all_data(initval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%I) First, collect an image series. TIP: run this once, then save the
%workspace including 'a3' and then comment out this section; saves lots of
%time with repaeated runs
    lbl3=strcat(initval.basepath,initval.driftfile);
    drift=dlmread(lbl3);
    lbl1=strcat(initval.basepath,initval.BFdatapath,initval.BFfiletemplate);
    aa = readtimeseries(lbl1,'tiff',[1 initval.maxfile]);              %Brightfielddata stack
    lbl2=strcat(initval.basepath,initval.FLdatapath,initval.FLfiletemplate);
    ff = readtimeseries(lbl2,'tiff',[1 initval.maxfile]);              %Fluorescencedata stack
    

