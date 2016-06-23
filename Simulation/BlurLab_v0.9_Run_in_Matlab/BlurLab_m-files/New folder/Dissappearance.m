function [X0,Y0,Z0,I0,L0,A0,R0,pts] = ...
    Dissappearance(X0,Y0,Z0,I0,L0,A0,R0,pts,CDis)

    Diss_cells = find(rand(pts,1) < CDis);

    X0(Diss_cells) = [];
    Y0(Diss_cells) = [];
    Z0(Diss_cells) = [];
    I0(Diss_cells) = [];
    L0(Diss_cells) = [];
    A0(Diss_cells) = [];
    R0(Diss_cells) = [];

    pts = length(X0);
end