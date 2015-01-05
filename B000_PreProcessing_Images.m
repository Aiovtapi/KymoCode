clc;
clear;

% for i = 23:30
    %strcat('Processing: Measurement ',' ',int2str(i),' of',int2str(10))
    
    %MeasurementNr = sprintf('%03i',i);
    ExpName = strcat('20131111_PAmCherry_DoubleIntensity_2mW/Slide1/Measurement_',MeasurementNr);
    initval=A001_Images_Set_Experiment(expno);
    
    strcat('RollingBalling')
    A000_RollingBall_FluorescenceImg(ExpName);
    strcat('Illuminiation correction')
    A000_IlluminationCorrection(ExpName);
%     strcat('PH and FL alignment')
%     %A000_PhaseContrast_Alignment(ExpName);
%end