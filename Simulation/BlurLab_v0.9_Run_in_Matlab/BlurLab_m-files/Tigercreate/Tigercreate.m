%Mark Shui Hu
%Random 3D Diffusion of particles

function Tigercreate(nframes,tif,pts,meanI,Lx,Ly,Lz,ini)
   
    %% Set initial values
    SigmaDperc = 0.10;
    
    % Define Growth of cell
    OLx = Lx / 2; 
    
    % Set initial number of blinking cells
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
    npts = pts + 1;
    J0 = I0;
    
    % Find Diffusion constant for each dimension
    [X, Y, D(1), D(2)] = TDDC;
    Dframes = round(X*(nframes - 1) + 1)';
    Dpart = Y'/100;
    
    % Get all values between the clicked points form TDDC
    for i = 1:(numel(Dframes)-1)
        Co = (Dpart(i+1) - Dpart(i)) / (Dframes(i+1) - Dframes(i));
        betx = ((Dframes(i) + 1) : (Dframes(i + 1) - 1))';
        cox = betx - Dframes(i);
        bety = Dpart(i) + cox * Co;
        
        Dframes = [Dframes; betx];
        Dpart = [Dpart; bety];
    end
    [~, sorti] = sort(Dframes);
    Dpart = Dpart(sorti);
    
    % Get the mean squared distance per frame for 2D Diffusion
    Dac = sqrt(4*D*tif); 
    SigmaD = Dac*SigmaDperc;
    
    % Select a part of spots to behave according the the first Diffusion
    % coefficient
    Allpts = (1:pts)';
    Dpts = ismember(Allpts,datasample(Allpts,round(Dpart(1)*pts),'Replace',false));
    
    % Get times travelled with current direction
    Framesec = 1/tif;
    T0 = Framesec*rand(pts,1);
    
    % Get initial change in direction
    % Get number of spots in each DC and get logicals
    d1pts = numel(find(Dpts));
    d2pts = pts - d1pts;
    Dpts = logical(Dpts);
    Dpts2 = ~Dpts;
    
    % Get Change in position for each spot
    XYC = zeros(pts,2);
    % Get distance travelled in both X and Y direction
    XYC(Dpts,:)  = 2*(randi([0,1],d1pts,2)-1/2) .* normrnd(Dac(1),SigmaD(1),d1pts,2);
    XYC(Dpts2,:) = 2*(randi([0,1],d2pts,2)-1/2) .* normrnd(Dac(2),SigmaD(2),d2pts,2);
    
    % Create empty output variables
    %   Xout: X position
    %   Yout: Y position
    %   Zout: Z position
    %   Iout: Intensity value
    %   Fout: Frame number
    %   Lout: Label; spot number
    
    [Xout,Yout,Zout,Iout,Fout,Lout] = deal([]);

    % Define growth of cell / box
    dLx = OLx / nframes;
    TLx = OLx;
    
    ini.Dac = Dac;
    ini.SigmaD = SigmaD;
    ini.Framesec = Framesec;
%%
    for i_fr = 1:nframes
        
        TLx = TLx + dLx;
        
        % Save outputs after each frame has passed for the output file
        NewL = [];
        Xout = [Xout; X0];
        Yout = [Yout; Y0];
        Zout = [Zout; Z0];
        Iout = [Iout; J0];
        Fout = [Fout; ones(pts,1)*i_fr];
        Lout = [Lout; L0];
        
        % Check whether the porportion of D1/D2 still holds
        % If not, some spots will be selected to change to the other
        % Diffusion constant
        if (sum(Dpts)/pts - Dpart(i_fr)) > 0.01
            % Too much in D1
            num = round((sum(Dpts)/pts - Dpart(i_fr))*pts);
            turn = datasample(find(Dpts),num);
            Dpts(turn) = false;
        elseif (Dpart(i_fr) - sum(Dpts)/pts) > 0.01
            % Too much in D2
            num = round((Dpart(i_fr) - sum(Dpts)/pts)*pts);
            turn = datasample(find(~Dpts),num);
            Dpts(turn) = true;
        else
            turn = [];
        end
        % find all spots that changed DC
        DCchange = false(pts,1);
        DCchange(turn) = true;
        
        % Find all spots travelled long enough for a change in direction
        % and reset counter
        Dchangeidx = T0 > Framesec;
        T0 = T0 + 1;
        T0(Dchangeidx) = T0(Dchangeidx) - Framesec;
        
        % Find spots with change in direction or DC, with both DC
        Allchange = DCchange | Dchangeidx;
        Changepts1 = Dpts & Allchange;
        Changepts2 = ~Dpts & Allchange;
        d1pts = sum(Changepts1);
        d2pts = sum(Changepts2);
        
        % Get new change in X and Y if neccessary
        XYC(Changepts1,:) = 2*(randi([0,1],d1pts,2)-1/2) .* normrnd(Dac(1),SigmaD(1),d1pts,2);
        XYC(Changepts2,:) = 2*(randi([0,1],d2pts,2)-1/2) .* normrnd(Dac(2),SigmaD(2),d2pts,2);
        
        % Apply change to positions
        X0 = X0 + XYC(:,1);
        Y0 = Y0 + XYC(:,2);
        
        % Find spots that have gone out of bounds
        Xstay = X0<0 | X0>TLx;
        Ystay = Y0<0 | Y0>Ly;
        
        % Get new change in X and Y
        XYC(Xstay,1) = -XYC(Xstay,1);
        XYC(Ystay,2) = -XYC(Ystay,2);
        
        % Limit Movement of out-of-bound spots
        X0(Xstay) = X0(Xstay) + 2*XYC(Xstay,1);
        Y0(Ystay) = Y0(Ystay) + 2*XYC(Ystay,2);
        
        % Change in Intensity
        Ichange = normrnd(1,ini.SigmaI,pts,1);
        I0 = I0.*Ichange;

        %% Merging
        
        % XYC,Dac,SigmaD T0, Framesec
        
        if ini.tog_merging == 1;
            [X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts,npts,NewL] = ...
                Tig_merging(X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts,npts,NewL,ini,i_fr);
        end

        %% Splitting
        if ini.tog_split == 1;
            [X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts,npts,NewL] = ...
                Tig_splitting(X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts,npts,NewL,ini,i_fr);
        end
        
        %% Spot Disappearence
        if ini.tog_spotdis == 1;
            [X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts] = ...
                Tig_Dissappearance(X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts,ini);
        end

        %% Spot Appearence
        if ini.tog_spotapp == 1; 
            [X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts,npts,NewL] = ...
                Tig_Appearance(X0,Y0,Z0,XYC,T0,I0,L0,Dpts,npts,meanI,NewL,ini,TLx,Ly,Lz);
        end

        %% Blinking
        if ini.tog_blinking == 1;

            % Remove gone spots
            D0(~ismember(D0(:,1),L0),:) = [];

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
            Ridx = rand(Nspots,1) < ini.Cblink2;
            D0(Ridx,:) = [];

            % sort blinking spots        
            [~, sortd] = sort(D0(:,1));
            D0 = D0(sortd,:);
        end
        
        %%
        
        % Remove low intensity spots      
        [X0,Y0,Z0,XYC,T0,I0,L0,Dpts,pts] = ...
            Tig_MinIntensity(X0,Y0,Z0,XYC,T0,I0,L0,Dpts,ini);
        
        % Compute real intensity (with blinking) (using histc to take into
        % account the appearence and dissapearence of spots
        J0 = I0;
        if numel(D0)>0
            [~,blinkspots] = histc(D0(:,1),L0);
            blinkI = I0(blinkspots).*D0(:,2);
            J0(blinkspots) = blinkI;        
        end
        
        clear NewL sortd 
    end
    
    %create file name
    tmp1=clock;
    fname=['BlurLab_rand_output_' date '_' num2str(tmp1(4)) '-' num2str(tmp1(5)) '-' num2str(tmp1(6))...
        '_pts-' num2str(pts) '_meanI-' num2str(meanI) '_D-' num2str(D(1)) '_Lx-' num2str(Lx) '_Ly-' ...
        num2str(Ly) '_Lz-' num2str(Lz) '.txt'];
    
    % write output file
    oldfolder = pwd;
    cd(ini.Tigfolder)
    BlurLab_text(Xout,Yout,Zout,Iout,Fout,Lout,fname)
    cd(oldfolder);
end



