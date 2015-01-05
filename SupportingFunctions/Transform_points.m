function point=Transform_points(point1D, initval, presets);
      %This function translates a coordinate along a 1D axis
      %to image-based coordinates. The 1D axis is defined by the original
      %channel selection, stored in 'presets' and 'initval' 
      
      r0=presets.twopoints(1,1); c0=presets.twopoints(1,2);
      r1=presets.twopoints(2,1); c1=presets.twopoints(2,2);
      fracpos=1-point1D/initval.kymolength;
      r2=r0+(r1-r0)*fracpos;
      c2=c0+(c1-c0)*fracpos;
      point=[r2 c2];  