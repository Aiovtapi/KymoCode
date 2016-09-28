
function [X0,Y0,Z0,I0,L0,D0,pts] = ...
    Tig_Dissappearance(X0,Y0,Z0,I0,L0,D0,pts,CDis)

    Diss_spots = find(rand(pts,1) < CDis);

    X0(Diss_spots) = [];
    Y0(Diss_spots) = [];
    Z0(Diss_spots) = [];
    I0(Diss_spots) = [];
    L0(Diss_spots) = [];
    D0(Diss_spots) = [];

    pts = length(X0);
end