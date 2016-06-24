%Mark Shui Hu
%Random 3D Diffusion of particles

function Tigercreate(nframes,pts,meanI,D,Lx,Ly,Lz,ini)
    %create file name
    tmp1=clock;
    fname=['BlurLab_rand_output_' date '_' num2str(tmp1(4)) '-' num2str(tmp1(5)) '-' num2str(tmp1(6))...
        '_pts-' num2str(pts) '_meanI-' num2str(meanI) '_D-' num2str(D) '_Lx-' num2str(Lx) '_Ly-' num2str(Ly) '_Lz-' num2str(Lz) '.txt'];


    %% Set initial values

    D0(:,1) = round(pts*rand(ini.Nblink,1));
    D0(:,2) = round(rand(ini.Nblink,1));
    D0(:,3) = round(ini.Tblink*rand(ini.Nblink,1));
    D0(:,4) = raylrnd(ini.Tblink,ini.Nblink,1);
    [~, sortD] = sort(D0(:,1));
    D0 = D0(sortD,:);

    X0 = Lx*rand(pts,1);
    Y0 = Ly*rand(pts,1);
    Z0 = Lz*rand(pts,1);
    I0 = poissrnd(meanI,pts,1);
    L0 = (1:pts)';
    A0 = 2*pi*rand(pts,1);
    R0 = normrnd(0,sqrt(2*D),pts,1);
    npts = pts + 1;
    J0 = I0;

    [Xout,Yout,Zout,Iout,Fout,Lout] = deal([]);

    %% 

    for i_fr = 1:nframes

        NewL = [];
        Xout = [Xout; X0];
        Yout = [Yout; Y0];
        Zout = [Zout; Z0];
        Iout = [Iout; J0];
        Fout = [Fout; ones(pts,1)*i_fr];
        Lout = [Lout; L0];

        X0 = X0 + R0.*sin(A0);
        Y0 = Y0 + R0.*cos(A0);

        [~,blinkspots] = histc(D0(:,1),L0);
        blinkI = I0(blinkspots).*D0(:,2);
        J0(blinkspots) = blinkI;

        %% Merging
        if ini.tog_merging == 1;
            [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
                Tig_merging(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,ini.MergeD,ini.CMerge,i_fr);
        end

        %% Splitting
        if ini.tog_split == 1;
            [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
                Tig_splitting(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,ini.CSplit,ini.DSplit,i_fr,meanI);

        end

        %% Spot Disappearence
        if ini.tog_spotdis == 1;
            [X0,Y0,Z0,I0,L0,A0,R0,pts] = ...
                Tig_Dissappearance(X0,Y0,Z0,I0,L0,A0,R0,pts,ini.CDis);
        end

        %% Spot Appearence
        if ini.tog_spotapp == 1; 
            [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
                Tig_Appearance(X0,Y0,Z0,I0,L0,A0,R0,npts,meanI,NewL,ini.Apois,Lx,Ly,Lz,D);
        end

        %% Blinking
        if ini.tog_blinking == 1;

            % Remove gone spots
            D0(find(ismember(D0(:,1),L0)==0),:) = [];

            % Add new blinking spots
            Nnspots = length(NewL);
            if ~Nnspots == 0;
                Newidx = find(rand(Nnspots,1) < ini.Cblink1);
                Dnew(:,1) = NewL(Newidx);
                Dnew(:,2) = round(rand(numel(Newidx),1));
                Dnew(:,3) = round(ini.Tblink*rand(numel(Newidx),1));
                Dnew(:,4) = raylrnd(ini.Tblink,numel(Newidx),1);

                D0 = [D0; Dnew];
            end

            % Add frame to blinkcounter
            D0(:,3) = D0(:,3) + 1;

            blinkaction = find(D0(:,3)>D0(:,4));
            ini.Nblink = numel(blinkaction);

            % Blinkaction
            for i_blink = blinkaction;
                if D0(i_blink,2) == 1;
                    D0(i_blink,2) = 0;
                elseif D0(i_blink,2) == 0;
                    D0(i_blink,2) = 1;
                end
                D0(i_blink,3) = 0;
                D0(i_blink,4) = raylrnd(ini.Tblink);
            end

            % Remove some blinking spots
            Nspots = size(D0,1);
            Ridx = find(rand(Nspots,1) < ini.Cblink2);
            D0(Ridx,:) = [];

            % sort blinking spots        
            [~, sortd] = sort(D0(:,1));
            D0 = D0(sortd,:);
        end

        J0 = I0;

        Rchange = normrnd(1,ini.SigmaR,pts,1);
        Achange = normrnd(1,ini.SigmaA,pts,1);
        Ichange = normrnd(1,ini.SigmaI,pts,1);

        R0 = R0.*Rchange;
        A0 = A0.*Achange;
        I0 = I0.*Ichange;


        clear NewL sortd 
    end

    %I know this isn't quite right
    Xout(Xout<0)=0;
    Xout(Xout>Lx)=Lx;
    Yout(Yout<0)=0;
    Yout(Yout>Ly)=Ly;
    Zout(Zout<0)=0;
    Zout(Zout>Lz)=Lz;

    % write output file
    BlurLab_text(Xout,Yout,Zout,Iout,Fout,Lout,fname)
end

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

    
function [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
    Tig_Appearance(X0,Y0,Z0,I0,L0,A0,R0,npts,meanI,NewL,Apois,Lx,Ly,Lz,D)

    N_app = round(poissrnd(Apois));

    L_new = (npts:npts+N_app-1)';

    X0 = [X0; Lx*rand(N_app,1)];
    Y0 = [Y0; Ly*rand(N_app,1)];
    Z0 = [Z0; Lz*rand(N_app,1)];
    I0 = [I0; poissrnd(meanI,N_app,1)];
    L0 = [L0; L_new];
    A0 = [A0; 2*pi*rand(N_app,1)];
    R0 = [R0; normrnd(0,sqrt(2*D),N_app,1)];
    npts = npts + N_app;

    pts = length(X0);
    NewL = [NewL; L_new];
end

function [X0,Y0,Z0,I0,L0,A0,R0,pts] = ...
    Tig_Dissappearance(X0,Y0,Z0,I0,L0,A0,R0,pts,CDis)

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
