function [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
    Tig_merging(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,MergeD,CMerge,i_fr)

    [Xa,Xb] = meshgrid(X0);
    [Ya,Yb] = meshgrid(Y0);
    Dist = sqrt(abs(Xa-Xb).^2+abs(Ya-Yb).^2);
    Smallidx = find(Dist<MergeD);
    Mergeidx = Smallidx(~ismember(Smallidx,find(diag(ones(1,pts)))));
    [Rowidx, Colidx] = ind2sub(size(Dist),Mergeidx);
    [~,Samep] = ismember(Colidx,Rowidx);
    duplicates = ones(size(Mergeidx));
    Dupidx = zeros(length(Mergeidx)/2,1);
    l = 1;

    for i_dup = 1:length(Mergeidx)
        if duplicates(i_dup) == 1;
            Dupidx(l) = Samep(i_dup);
            l = l + 1;
            duplicates(Samep(i_dup)) = 0;
        end
    end  

    Colidx(Dupidx) = [];
    Rowidx(Dupidx) = [];
    Cont_merge = find(rand(length(Colidx),1) < CMerge);
    Merge_cells = [Rowidx(Cont_merge), Colidx(Cont_merge)];

    N_merge = size(Merge_cells,1);
    [X_new,Y_new,Z_new,I_new,F_new,L_new,R_new,A_new] = deal(zeros(N_merge,1));

    for i_mrg = 1:N_merge;
        cell_a = Merge_cells(i_mrg,1);
        cell_b = Merge_cells(i_mrg,2);

        X_new(i_mrg) = (X0(cell_a)+X0(cell_b))/2;
        Y_new(i_mrg) = (Y0(cell_a)+Y0(cell_b))/2;
        Z_new(i_mrg) = (Z0(cell_a)+Z0(cell_b))/2;
        I_new(i_mrg) = I0(cell_a)+I0(cell_b);
        F_new(i_mrg) = i_fr;
        L_new(i_mrg) = npts;

        Mx_a = R0(cell_a)*sin(A0(cell_a))*I0(cell_a);
        Mx_b = R0(cell_b)*sin(A0(cell_b))*I0(cell_b);
        My_a = R0(cell_a)*cos(A0(cell_a))*I0(cell_a);
        My_b = R0(cell_b)*cos(A0(cell_b))*I0(cell_b);

        Rx_new = (Mx_a + Mx_b)/I_new(i_mrg);
        Ry_new = (My_a + My_b)/I_new(i_mrg);
        A_new(i_mrg) = tan(Rx_new/Ry_new);
        R_new(i_mrg) = sqrt(Rx_new^2 + Ry_new^2); 

        npts = npts + 1;
    end

    Rcells = Merge_cells(:);
    X0(Rcells) = [];    X0 = [X0; X_new];
    Y0(Rcells) = [];    Y0 = [Y0; Y_new];
    Z0(Rcells) = [];    Z0 = [Z0; Z_new];
    I0(Rcells) = [];    I0 = [I0; I_new];
    L0(Rcells) = [];    L0 = [L0; L_new];
    A0(Rcells) = [];    A0 = [A0; A_new];
    R0(Rcells) = [];    R0 = [R0; R_new];

    pts = length(X0);
    NewL = [NewL; L_new];
end