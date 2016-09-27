
function A0 = Tig_wallcol(Xstay,Ystay,A0)
    
    XYstay = Xstay & Ystay;
    Xstay = logical(Xstay - XYstay);
    Ystay = logical(Ystay - XYstay); 
    
    XA = A0(Xstay);
    YA = A0(Ystay);
    XYA = A0(XYstay);
    
    for i = 1:numel(XA);
        ag = XA(i);
        if ag <= pi
            ag = pi - ag; 
        else
            ag = 3*pi - ag; 
        end
        XA(i) = ag;         
    end
    
    YA = 2*pi - YA; 
    XYA = XYA - pi;
    
    A0(Xstay) = XA;
    A0(Ystay) = YA;
    A0(XYstay) = XYA;
end

