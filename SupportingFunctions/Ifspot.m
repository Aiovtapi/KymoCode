 function div=Ifspot(Divsoi,tr);
 %This function detrmines if a short subsection contains a spot
     div=0;
     
     if range(Divsoi)>tr,
         div=1;,
     end
  
  end