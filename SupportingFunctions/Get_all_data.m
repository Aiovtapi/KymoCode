  function [aa,ff,drift]=Get_all_data(initval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%I) First, collect an image series. TIP: run this once, then save the
%workspace including 'a3' and then comment out this section; saves lots of
%time with repeated runs

    lbl3=strcat(initval.basepath,initval.driftfile);
    
    if exist(lbl3,'file')==0
        fileID=fopen(lbl3,'a');
        fprintf(fileID,'0');
        fclose(fileID);
        drift=dlmread(lbl3);
    elseif exist(lbl3,'file')==2
    drift=dlmread(lbl3);
    end
    
    lbl1=strcat(initval.basepath,initval.BFdatapath,initval.BFfiletemplate);
    aa = readtimeseries(lbl1,'tiff',[1 initval.maxfile]);              %Brightfielddata stack
    lbl2=strcat(initval.basepath,initval.FLdatapath,initval.FLfiletemplate);
    ff = readtimeseries(lbl2,'tiff',[1 initval.maxfile]);              %Fluorescencedata stack
    

