
function [X0,Y0,Z0,I0,L0,A0,R0,D0,pts] = ...
    Tig_MinIntensity(X0,Y0,Z0,I0,L0,A0,R0,D0,MinI)

    Diss_spots = find(I0 < MinI);

    X0(Diss_spots) = [];
    Y0(Diss_spots) = [];
    Z0(Diss_spots) = [];
    I0(Diss_spots) = [];
    L0(Diss_spots) = [];
    A0(Diss_spots) = [];
    R0(Diss_spots) = [];
    D0(Diss_spots) = [];

    pts = length(X0);
end
