function [X0,Y0,Z0,XYC,T0,I0,L0,D0,pts,npts,NewL] = ...
    Tig_splitting(X0,Y0,Z0,XYC,T0,I0,L0,D0,pts,npts,NewL,ini,i_fr)

    % Get random splitting spots
    Split_spots = find(rand(pts,1) < ini.CSplit);

    % Get number of splitting spots
    N_split = size(Split_spots,1);
    [X_new,Y_new,Z_new,I_new,F_new,L_new,T_new] = deal(zeros(N_split*2,1));
    D_new = false(N_split*2,1);
    XYC_new = zeros(N_split,2);

    for i_splt = 1:N_split;
        % Get spot number
        spot_old = Split_spots(i_splt);
        
        % Get new intensity values
        I_old = I0(spot_old);            
        if ini.tog_dimer == 0
            I_a = poissrnd(I_old/2);
        else
            I_a = I_old/2;
        end
        I_b = I_old - I_a;
        
        % Save positions of new spots
        X_new(2*i_splt-1:2*i_splt) = X0(spot_old);
        Y_new(2*i_splt-1:2*i_splt) = Y0(spot_old);
        Z_new(2*i_splt-1:2*i_splt) = Z0(spot_old);
        
        % Save other properties of new spots
        I_new(2*i_splt-1) = I_a;
        I_new(2*i_splt) = I_b;
        F_new(2*i_splt-1:2*i_splt) = i_fr;
        L_new(2*i_splt-1) = npts;
        L_new(2*i_splt) = npts + 1;
        D_new(2*i_splt-1:2*i_splt) = logical(D0(spot_old));
        T_new(2*i_splt-1:2*i_splt) = T0(spot_old);
        
        % Set new direction
        if D0(spot_old)
            Daci = 1;
        else
            Daci = 2;
        end
        XYC_new(2*i_splt-1:2*i_splt,:) = 2*(randi([0,1],2,2)-1/2) .* ...
            normrnd(ini.Dac(Daci),ini.SigmaD(Daci),2,2);

        npts = npts + 2;
    end
    
    % Save values to matrices
    X0(Split_spots) = [];    X0 = [X0; X_new];
    Y0(Split_spots) = [];    Y0 = [Y0; Y_new];
    Z0(Split_spots) = [];    Z0 = [Z0; Z_new];
    I0(Split_spots) = [];    I0 = [I0; I_new];
    L0(Split_spots) = [];    L0 = [L0; L_new];
    D0(Split_spots) = [];    D0 = [D0; D_new];
    T0(Split_spots) = [];    T0 = [T0; T_new];
    XYC(Split_spots,:) = [];   XYC = [XYC; XYC_new];
    
    pts = length(X0);
    NewL = [NewL; L_new];
end

