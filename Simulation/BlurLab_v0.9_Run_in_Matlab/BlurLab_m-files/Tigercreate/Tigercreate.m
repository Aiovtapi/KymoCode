%Mark Shui Hu
%Random 3D Diffusion of particles

function Tigercreate(nframes,pts,meanI,D,Lx,Ly,Lz,ini)
    %create file name
    tmp1=clock;
    fname=['BlurLab_rand_output_' date '_' num2str(tmp1(4)) '-' num2str(tmp1(5)) '-' num2str(tmp1(6))...
        '_pts-' num2str(pts) '_meanI-' num2str(meanI) '_D-' num2str(D) '_Lx-' num2str(Lx) '_Ly-' num2str(Ly) '_Lz-' num2str(Lz) '.txt'];


    %% Set initial values
    SigmaDperc = 0.10;
    
    
    OLx = Lx / 2; 
    Nblink2 = round(pts*ini.Nblink/100);

    % D0 Matrix contains only the blinking spots with in...
    %   colum 1: The index of the spot.
    %   colum 2: The current state of the spot (on/off)
    %   colum 3: The time spent the this state
    %   colum 4: The target time the spot will stay in its current state
    
    D0(:,1) = ceil(pts*rand(Nblink2,1));
    D0(:,2) = round(rand(Nblink2,1));
    D0(:,3) = round(ini.Tblink*rand(Nblink2,1));
    D0(:,4) = raylrnd(ini.Tblink,Nblink2,1);
    
    % Sort D0 to the indices of the spots
    [~, sortD] = sort(D0(:,1));
    D0 = D0(sortD,:);
    
    % Create initial values of position, intensity and movement
    X0 = OLx*rand(pts,1);
    Y0 = Ly*rand(pts,1);
    Z0 = Lz*rand(pts,1);
    
    I0 = poissrnd(meanI,pts,1);
    L0 = (1:pts)';
    A0 = 2*pi*rand(pts,1);
    R0 = normrnd(0,sqrt(2*D),pts,1);
    npts = pts + 1;
    J0 = I0;
    
    Dim = 0;
    if Lx > 1; Dim = Dim + 1; end
    if Ly > 1; Dim = Dim + 1; end
    if Lz > 1; Dim = Dim + 1; end
    Dac = sqrt(D^2/Dim);
    SigmaD = Dac/SigmaDperc;
    % Create empty output variables
    %   Xout: X position
    %   Yout: Y position
    %   Zout: Z position
    %   Iout: Intensity value
    %   Fout: Frame number
    %   Lout: Label; spot number
    
    [Xout,Yout,Zout,Iout,Fout,Lout] = deal([]);

    %% 
    
    dLx = OLx / nframes;
    TLx = OLx;

    for i_fr = 1:nframes
        
        TLx = TLx + dLx;

        NewL = [];
        Xout = [Xout; X0];
        Yout = [Yout; Y0];
        Zout = [Zout; Z0];
        Iout = [Iout; J0];
        Fout = [Fout; ones(pts,1)*i_fr];
        Lout = [Lout; L0];

%         % Compute new XY position
%         X0 = X0 + R0.*sin(A0);
%         Y0 = Y0 + R0.*cos(A0);
%         
%         % Find out-of-bound spots
%         Xstay = X0<0 | X0>TLx;
%         Ystay = Y0<0 | Y0>Ly;
%         
%         % Limit Movement of out-of-bound spots
%         X0(Xstay) = X0(Xstay) - R0(Xstay).*sin(A0(Xstay));
%         Y0(Ystay) = Y0(Ystay) - R0(Ystay).*cos(A0(Ystay));
%         
%         % Calculate new angles after wall collision for out-of-bound spots
%         A0 = wrapTo2Pi(A0);
%         A0 = Tig_wallcol(Xstay,Ystay,A0);
%         Allstay = Xstay | Ystay;
%         
%         % Compute new YX position after wall collision
%         X0(Allstay) = X0(Allstay) + R0(Allstay).*sin(A0(Allstay));
%         Y0(Allstay) = Y0(Allstay) + R0(Allstay).*cos(A0(Allstay));

        
        Xchange = normrnd(Dac,SigmaD,pts,1);
        Ychange = normrnd(Dac,SigmaD,pts,1);
        Zchange = normrnd(Dac,SigmaD,pts,1);

        X0 = X0 + Xchange;
        Y0 = Y0 + Ychange;
        Z0 = Z0 + Zchange;
        
        Xstay = X0<0 | X0>TLx;
        Ystay = Y0<0 | Y0>Ly;
        Zstay = Z0<0 | Z0>Lz;
        
        % Limit Movement of out-of-bound spots
        X0(Xstay) = X0(Xstay) - Xchange(Xstay);
        Y0(Ystay) = Y0(Ystay) - Ychange(Ystay);
        Z0(Zstay) = Z0(Zstay) - Zchange(Zstay);
        
        % Compute real intensity (with blinking) (using histc to take into
        % account the appearence and dissapearence of spots
        if numel(D0)>0
            [~,blinkspots] = histc(D0(:,1),L0);
            blinkI = I0(blinkspots).*D0(:,2);
            J0(blinkspots) = blinkI;        
        end


        %% Merging
        if ini.tog_merging == 1;
            [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
                Tig_merging(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,ini.MergeD,ini.CMerge,i_fr);
        end

        %% Splitting
        if ini.tog_split == 1;
            [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
                Tig_splitting(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,ini.CSplit,ini.DSplit,ini.tog_dimer,i_fr,meanI);
        end


        %% Spot Disappearence
        if ini.tog_spotdis == 1;
            [X0,Y0,Z0,I0,L0,A0,R0,pts] = ...
                Tig_Dissappearance(X0,Y0,Z0,I0,L0,A0,R0,pts,ini.CDis);
        end


        %% Spot Appearence
        if ini.tog_spotapp == 1; 
            [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
                Tig_Appearance(X0,Y0,Z0,I0,L0,A0,R0,npts,meanI,NewL,ini.Apois,TLx,Ly,Lz,D);
        end

        %% Blinking
        if ini.tog_blinking == 1;

            % Remove gone spots
            D0(find(ismember(D0(:,1),L0)==0),:) = [];

            % Add new blinking spots
            Nnspots = length(NewL);
            if ~Nnspots == 0;
                Newidx = find(rand(Nnspots,1) < ini.Cblink1);
                if numel(Newidx)>0
                    
                    Dnew(:,1) = NewL(Newidx);
                    Dnew(:,2) = round(rand(numel(Newidx),1));
                    Dnew(:,3) = round(ini.Tblink*rand(numel(Newidx),1));
                    Dnew(:,4) = raylrnd(ini.Tblink,numel(Newidx),1);

                    D0 = [D0; Dnew];
                    clear Dnew
                end
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

        Rchange = normrnd(1,ini.SigmaR,pts,1);
        Achange = normrnd(1,ini.SigmaA,pts,1);
        Ichange = normrnd(1,ini.SigmaI,pts,1);

        R0 = R0.*Rchange;
        A0 = A0.*Achange;
        I0 = I0.*Ichange;
        
%       Remove low intensity spots      
        [X0,Y0,Z0,I0,L0,A0,R0,pts] = ...
            Tig_MinIntensity(X0,Y0,Z0,I0,L0,A0,R0,ini.MinI);
        J0 = I0;

        clear NewL sortd 
    end

    % write output file
    oldfolder = pwd;
    cd(ini.Tigfolder)
    BlurLab_text(Xout,Yout,Zout,Iout,Fout,Lout,fname)
    cd(oldfolder);
end

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
    Merge_spots = [Rowidx(Cont_merge), Colidx(Cont_merge)];

    N_merge = size(Merge_spots,1);
    [X_new,Y_new,Z_new,I_new,F_new,L_new,R_new,A_new] = deal(zeros(N_merge,1));

    for i_mrg = 1:N_merge;
        cell_a = Merge_spots(i_mrg,1);
        cell_b = Merge_spots(i_mrg,2);

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

    Rspots = Merge_spots(:);
    X0(Rspots) = [];    X0 = [X0; X_new];
    Y0(Rspots) = [];    Y0 = [Y0; Y_new];
    Z0(Rspots) = [];    Z0 = [Z0; Z_new];
    I0(Rspots) = [];    I0 = [I0; I_new];
    L0(Rspots) = [];    L0 = [L0; L_new];
    A0(Rspots) = [];    A0 = [A0; A_new];
    R0(Rspots) = [];    R0 = [R0; R_new];

    pts = length(X0);
    NewL = [NewL; L_new];
end

function [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
    Tig_splitting(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,CSplit,DSplit,tog_dimer,i_fr,meanI)

    Split_spots = find(rand(pts,1) < CSplit);

    N_split = size(Split_spots,1);
    [X_new,Y_new,Z_new,I_new,F_new,L_new,R_new,A_new] = deal(zeros(N_split*2,1));

    for i_splt = 1:N_split;
        cell_old = Split_spots(i_splt);

        I_old = I0(cell_old);            
        Mx_old = R0(cell_old)*sin(A0(cell_old))*I0(cell_old);
        My_old = R0(cell_old)*cos(A0(cell_old))*I0(cell_old);

        if tog_dimer == 0
            Mx_a = normrnd(Mx_old/2,sqrt(2*DSplit)*sqrt(meanI));
            My_a = normrnd(My_old/2,sqrt(2*DSplit)*sqrt(meanI));
            I_a = poissrnd(I_old/2);
        else
            Mx_a = Mx_old/2;
            My_a = My_old/2;
            I_a = I_old/2;
        end
        
        Mx_b = Mx_old - Mx_a;
        My_b = My_old - My_a;
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

    X0(Split_spots) = [];    X0 = [X0; X_new];
    Y0(Split_spots) = [];    Y0 = [Y0; Y_new];
    Z0(Split_spots) = [];    Z0 = [Z0; Z_new];
    I0(Split_spots) = [];    I0 = [I0; I_new];
    L0(Split_spots) = [];    L0 = [L0; L_new];
    A0(Split_spots) = [];    A0 = [A0; A_new];
    R0(Split_spots) = [];    R0 = [R0; R_new];

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

    Diss_spots = find(rand(pts,1) < CDis);

    X0(Diss_spots) = [];
    Y0(Diss_spots) = [];
    Z0(Diss_spots) = [];
    I0(Diss_spots) = [];
    L0(Diss_spots) = [];
    A0(Diss_spots) = [];
    R0(Diss_spots) = [];

    pts = length(X0);
end

function [X0,Y0,Z0,I0,L0,A0,R0,pts] = ...
    Tig_MinIntensity(X0,Y0,Z0,I0,L0,A0,R0,MinI)

    Diss_spots = find(I0 < MinI);

    X0(Diss_spots) = [];
    Y0(Diss_spots) = [];
    Z0(Diss_spots) = [];
    I0(Diss_spots) = [];
    L0(Diss_spots) = [];
    A0(Diss_spots) = [];
    R0(Diss_spots) = [];

    pts = length(X0);
end
