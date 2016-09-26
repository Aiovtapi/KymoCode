function [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
    Tig_splitting(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,CSplit,DSplit,i_fr,meanI)

    Split_cells = find(rand(pts,1) < CSplit);

    N_split = size(Split_cells,1);
    [X_new,Y_new,Z_new,I_new,F_new,L_new,R_new,A_new] = deal(zeros(N_split*2,1));

    for i_splt = 1:N_split;
        cell_old = Split_cells(i_splt);

        I_old = I0(cell_old);            
        Mx_old = R0(cell_old)*sin(A0(cell_old))*I0(cell_old);
        My_old = R0(cell_old)*cos(A0(cell_old))*I0(cell_old);

        Mx_a = normrnd(Mx_old,sqrt(2*DSplit)*meanI);
        My_a = normrnd(My_old,sqrt(2*DSplit)*meanI);
        Mx_b = Mx_old - Mx_a;
        My_b = My_old - My_a;
        I_a = poissrnd(I_old);
        I_b = I_old - I_a;

        Rx_a = Mx_a/I_a;
        Ry_a = My_a/I_a;
        Rx_b = Mx_b/I_b;
        Ry_b = My_b/I_b;
        R_a = sqrt(Rx_a^2+Ry_a^2);
        R_b = sqrt(Rx_b^2+Ry_b^2);
        A_a = tan(Rx_a/Ry_a);
        A_b = tan(Rx_b/Ry_b);

        X_new(2*i_splt-1) = X0(cell_old) + Rx_a;
        X_new(2*i_splt) = X0(cell_old) + Rx_b;
        Y_new(2*i_splt-1) = Y0(cell_old) + Ry_a;
        Y_new(2*i_splt) = Y0(cell_old) + Ry_b;
        Z_new(2*i_splt-1:2*i_splt) = Z0(cell_old);

        I_new(2*i_splt-1) = I_a;
        I_new(2*i_splt) = I_b;
        F_new(2*i_splt-1:2*i_splt) = i_fr;
        L_new(2*i_splt-1) = npts;
        L_new(2*i_splt) = npts + 1;

        R_new(2*i_splt-1) = R_a;
        R_new(2*i_splt) = R_b;
        A_new(2*i_splt-1) = A_a;
        A_new(2*i_splt) = A_b;

        npts = npts + 2;
    end

    X0(Split_cells) = [];    X0 = [X0; X_new];
    Y0(Split_cells) = [];    Y0 = [Y0; Y_new];
    Z0(Split_cells) = [];    Z0 = [Z0; Z_new];
    I0(Split_cells) = [];    I0 = [I0; I_new];
    L0(Split_cells) = [];    L0 = [L0; L_new];
    A0(Split_cells) = [];    A0 = [A0; A_new];
    R0(Split_cells) = [];    R0 = [R0; R_new];

    pts = length(X0);
    NewL = [NewL; L_new];
end