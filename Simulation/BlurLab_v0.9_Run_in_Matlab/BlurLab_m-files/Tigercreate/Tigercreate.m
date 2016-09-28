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
    Dframes = round(X*(nframes - 1) + 1);
    Dpart = Y/100;
    
    % Get all values between the clicked points form TDDC
    for i = 1:(numel(Dframes)-1)
        Co = (Dpart(i+1) - Dpart(i)) / (Dframes(i+1) - Dframes(i));
        betx = (Dframes(i) + 1) : (Dframes(i + 1) - 1);
        cox = betx - Dframes(i);
        bety = Dpart(i) + cox * Co;
        
        Dframes = [Dframes, betx];
        Dpart = [Dpart, bety];
    end
    [~, sorti] = sort(Dframes);
    Dpart = Dpart(sorti);
    
    % Get the mean squared distance per frame for 2D Diffusion
    Dac = sqrt(4*D*tif); 
    SigmaD = Dac*SigmaDperc;
    
    % Select a part of spots to behave according the the first Diffusion
    % coefficient
    Allpts = 1:pts;
    Dpts = ismember(Allpts,datasample(Allpts,round(Dpart(1)*pts),'Replace',false));
    
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
    
%%
    for i_fr = 1:nframes

        clear Xchange Ychange 
        
        TLx = TLx + dLx;
        
        % Save outputs after each frame has passed for the output file
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
        end
        
        % Get number of spots in each DC and get logicals
        d1pts = numel(find(Dpts));
        d2pts = pts - d1pts;
        Dpts = logical(Dpts);
        Dpts2 = ~Dpts;
        
        % Get Change in position for each spot
        [Xchange,Ychange] = deal(zeros(pts,1));
        % Get direction in X and Y, either 1 or -1
        Xdir = 2*(randi([0,1],pts,1)-1/2);
        Ydir = 2*(randi([0,1],pts,1)-1/2);
        % Get distance travelled in both X and Y direction
        Xchange(Dpts) = normrnd(Dac(1),SigmaD(1),d1pts,1);
        Ychange(Dpts) = normrnd(Dac(1),SigmaD(1),d1pts,1);
        Xchange(Dpts2) = normrnd(Dac(2),SigmaD(2),d2pts,1);
        Ychange(Dpts2) = normrnd(Dac(2),SigmaD(2),d2pts,1);
        % Find final change in XY and apply them to initial positions
        Xchange = Xdir.*Xchange;
        Ychange = Ydir.*Ychange;
        X0 = X0 + Xchange;
        Y0 = Y0 + Ychange;
        
        % Find spots that have gone out of bounds
        Xstay = X0<0 | X0>TLx;
        Ystay = Y0<0 | Y0>Ly;
        
        % Get new change in X and Y
        Xchange(Xstay) = -Xchange(Xstay);
        Ychange(Ystay) = -Ychange(Ystay);
        
        % Limit Movement of out-of-bound spots
        X0(Xstay) = X0(Xstay) + 2*Xchange(Xstay);
        Y0(Ystay) = Y0(Ystay) + 2*Ychange(Ystay);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate R and A
        R0 = sqrt(Xchange.^2 + Ychange.^2);
        A0 = atan2(Ychange,Xchange);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
        % Compute real intensity (with blinking) (using histc to take into
        % account the appearence and dissapearence of spots
        if numel(D0)>0
            [~,blinkspots] = histc(D0(:,1),L0);
            blinkI = I0(blinkspots).*D0(:,2);
            J0(blinkspots) = blinkI;        
        end


        %% Merging
        
        if ini.tog_merging == 1;
            [X0,Y0,Z0,I0,L0,Dpts,pts,npts,NewL] = ...
                Tig_merging(X0,Y0,Z0,I0,L0,Dpts,pts,npts,NewL,ini.MergeD,ini.CMerge,i_fr);
        end

        %% Splitting
        if ini.tog_split == 1;
            [X0,Y0,Z0,I0,L0,A0,R0,Dpts,pts,npts,NewL] = ...
                Tig_splitting(X0,Y0,Z0,I0,L0,A0,R0,Dpts,pts,npts,NewL,ini.CSplit,ini.DSplit,ini.tog_dimer,i_fr,meanI);
        end
        

        %% Spot Disappearence
        if ini.tog_spotdis == 1;
            [X0,Y0,Z0,I0,L0,Dpts,pts] = ...
                Tig_Dissappearance(X0,Y0,Z0,I0,L0,Dpts,pts,ini.CDis);
        end

        %% Spot Appearence
        if ini.tog_spotapp == 1; 
            [X0,Y0,Z0,I0,L0,Dpts,pts,npts,NewL] = ...
                Tig_Appearance(X0,Y0,Z0,I0,L0,Dpts,npts,meanI,NewL,ini.Apois,TLx,Ly,Lz);
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

%         Rchange = normrnd(1,ini.SigmaR,pts,1);
%         Achange = normrnd(1,ini.SigmaA,pts,1);
        Ichange = normrnd(1,ini.SigmaI,pts,1);

%         R0 = R0.*Rchange;
%         A0 = A0.*Achange;
        I0 = I0.*Ichange;
        
%       Remove low intensity spots      
        [X0,Y0,Z0,I0,L0,Dpts,pts] = ...
            Tig_MinIntensity(X0,Y0,Z0,I0,L0,Dpts,ini.MinI);
        J0 = I0;

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



