
function [X0,Y0,Z0,XYC,T0,I0,L0,D0,pts] = ...
    Tig_MinIntensity(X0,Y0,Z0,XYC,T0,I0,L0,D0,ini)
    
    % Find spots that have intensities lower than set threshold
    Diss_spots = find(I0 < ini.MinI);

    X0(Diss_spots) = [];
    Y0(Diss_spots) = [];
    Z0(Diss_spots) = [];
    I0(Diss_spots) = [];
    L0(Diss_spots) = [];
    D0(Diss_spots) = [];
    T0(Diss_spots) = [];
    XYC(Diss_spots,:) = [];

    pts = length(X0);
end
