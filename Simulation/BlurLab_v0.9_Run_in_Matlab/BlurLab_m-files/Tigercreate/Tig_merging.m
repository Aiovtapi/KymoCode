function [X0,Y0,Z0,XYC,T0,I0,L0,D0,pts,npts,NewL] = ...
    Tig_merging(X0,Y0,Z0,XYC,T0,I0,L0,D0,pts,npts,NewL,ini,i_fr)
    
    % Create meshgrids to easily find distances between points
    [Xa,Xb] = meshgrid(X0);
    [Ya,Yb] = meshgrid(Y0);
    Dist = sqrt(abs(Xa-Xb).^2+abs(Ya-Yb).^2);
    
    % Find spots that are close to each other
    Smallidx = find(Dist<ini.MergeD);
    Mergeidx = Smallidx(~ismember(Smallidx,find(diag(ones(1,pts)))));
    [Rowidx, Colidx] = ind2sub(size(Dist),Mergeidx);
    [~,Samep] = ismember(Colidx,Rowidx);
    duplicates = ones(size(Mergeidx));
    Dupidx = zeros(length(Mergeidx)/2,1);
    l = 1;

    % Find duplicates and delete them
    for i_dup = 1:length(Mergeidx)
        if duplicates(i_dup) == 1;
            Dupidx(l) = Samep(i_dup);
            l = l + 1;
            duplicates(Samep(i_dup)) = 0;
        end
    end  
    Colidx(Dupidx) = [];
    Rowidx(Dupidx) = [];
    
    % Apply chance to merge event. Merge_spots contains the final list of
    % spots that will be merged
    Cont_merge = find(rand(length(Colidx),1) < ini.CMerge);
    Merge_spots = [Rowidx(Cont_merge), Colidx(Cont_merge)];
    
    % Get number of merging events and create variables
    N_merge = size(Merge_spots,1);
    [X_new,Y_new,Z_new,I_new,F_new,L_new,T_new] = deal(zeros(N_merge,1));
    XYC_new = zeros(N_merge,2);
    D_new = false(N_merge,1);

    for i_mrg = 1:N_merge;
        cell_a = Merge_spots(i_mrg,1);
        cell_b = Merge_spots(i_mrg,2);
        
        % Get new positions and define frame and spot number
        X_new(i_mrg) = (X0(cell_a)+X0(cell_b))/2;
        Y_new(i_mrg) = (Y0(cell_a)+Y0(cell_b))/2;
        Z_new(i_mrg) = (Z0(cell_a)+Z0(cell_b))/2;
        I_new(i_mrg) = I0(cell_a)+I0(cell_b);
        F_new(i_mrg) = i_fr;
        L_new(i_mrg) = npts;
        T_new(i_mrg) = (T0(cell_a)+T0(cell_b))/2;
        
        % Select the diffusion constant 
        if ~(D0(cell_a) == D0(cell_b))
            D_new(i_mrg) = randi([1,ini.numDC]);
        else
            D_new(i_mrg) = D0(cell_a);
        end
        
        % Set new direction change
        if D_new(i_mrg)
            Daci = 1;
        else
            Daci = 2;
        end
        XYC_new(i_mrg,:) = 2*(randi([0,1],1,2)-1/2) .* ...
            normrnd(ini.Dac(Daci),ini.SigmaD(Daci),1,2);

        npts = npts + 1;
    end
    
    % Delete old spots and add new spots
    Rspots = Merge_spots(:);
    X0(Rspots) = [];    X0 = [X0; X_new];
    Y0(Rspots) = [];    Y0 = [Y0; Y_new];
    Z0(Rspots) = [];    Z0 = [Z0; Z_new];
    I0(Rspots) = [];    I0 = [I0; I_new];
    L0(Rspots) = [];    L0 = [L0; L_new];
    D0(Rspots) = [];    D0 = [D0; D_new];
    T0(Rspots) = [];    T0 = [T0; T_new];
    XYC(Rspots,:) = [];   XYC = [XYC; XYC_new];

    pts = length(X0);
    NewL = [NewL; L_new];
end