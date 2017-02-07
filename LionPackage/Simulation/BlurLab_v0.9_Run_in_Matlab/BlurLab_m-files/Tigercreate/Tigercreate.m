%Mark Shui Hu
%Random 3D Diffusion of particles

function Tigercreate(nframes,tif,pts,meanI,Lx,Ly,Lz,ini)
   
    %% Set initial values
    SigmaDperc = 0.01;
    
    % Define Growth of cell
    if ini.Cellgrowth
        OLx = Lx / 2;
    else
        OLx = Lx;
    end
    
    % Set initial number of blinking cells
    Nblink2 = round(pts*ini.Nblink/100);

    % D0 Matrix contains only the blinking spots with in...
    %   colum 1: The index of the spot.
    %   colum 2: The current state of the spot (on/off)
    %   colum 3: The time spent the this state
    %   colum 4: The target time the spot will stay in its current state
    
    D0(:,1) = datasample(1:pts,Nblink2,'Replace',false)';
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
    
    if ini.tog_blinking == 1;
        if ~isempty(D0)
            [~,blinkspots] = histc(D0(:,1),L0);
            blinkI = I0(blinkspots).*D0(:,2);
            J0(blinkspots) = blinkI;        
        end
    end
    
    % Find Diffusion constant for each dimension
    numDC = ini.numDC;
    DC = ini.DC;
    DCY = ini.DCY;
    
    % Get the mean squared distance per frame for 2D Diffusion
    Dac = sqrt(4*DC*tif); 
    SigmaD = Dac*SigmaDperc;
    
    % Select a part of spots to behave according the the first Diffusion
    % Get initial change in direction
    % Get number of spots in each DC and get logicals
    Allpts = (1:pts)';
    
    d1pts = pts;            d2pts = 0;              d3pts = 0;
    Dpts1 = true(pts,1);    Dpts2 = false(pts,1);   Dpts3 = false(pts,1);
    
    if numDC > 1
        Dpts2 = ismember(Allpts,datasample(Allpts,round(DCY(1,2)*pts),'Replace',false));
        d2pts = sum(Dpts2);
        Dpts1 = xor(Dpts1,Dpts2);
        d1pts = sum(Dpts1);
    end
    if numDC > 2
        Partpts = Allpts(~Dpts2);
        ppts = numel(Partpts);
        if ~isempty(Partpts) && ~(DCY(1,2) == 1)
            Dpts3 = ismember(Allpts,datasample(Partpts,round(DCY(1,3)/(1-DCY(1,2))*ppts),'Replace',false));
            d3pts = sum(Dpts3);
            Dpts1 = xor(Dpts1,Dpts3);
            d1pts = sum(Dpts1);
        end
    end
    
    DCpts = double(Dpts1)*1 + double(Dpts2)*2 + double(Dpts3)*3;

    % Get times travelled with current direction
    Framesec = 1/tif;
    T0 = Framesec*rand(pts,1);
    
    % Get Change in position for each spot
    XYC = zeros(pts,2);
    % Get distance travelled in both X and Y direction
    XYC(Dpts1,:) = 2*(randi([0,1],d1pts,2)-1/2) .* normrnd(Dac(1),SigmaD(1),d1pts,2);
    XYC(Dpts2,:) = 2*(randi([0,1],d2pts,2)-1/2) .* normrnd(Dac(2),SigmaD(2),d2pts,2);
    XYC(Dpts3,:) = 2*(randi([0,1],d3pts,2)-1/2) .* normrnd(Dac(3),SigmaD(3),d3pts,2);
    
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
        if ini.Cellgrowth
            TLx = TLx + dLx;
        end
        
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
        DC2much = [(DCY(i_fr,1) - sum(DCpts == 1)/pts) > 0.01,...
                   (DCY(i_fr,2) - sum(DCpts == 2)/pts) > 0.01,...
                   (DCY(i_fr,3) - sum(DCpts == 3)/pts) > 0.01];
        
        totalturn = [];
        
        if DC2much(1) % Too little with DC1
            if      ~DC2much(2) && DC2much(3)
                Changeops = find(DCpts == 2);
            elseif  DC2much(2) && ~DC2much(3)
                Changeops = find(DCpts == 3);
            elseif  ~DC2much(2) && ~DC2much(3)
                Changeops = find(DCpts == 2 | DCpts == 3);
            end
            Changenum = round(DCY(i_fr,1)*pts - sum(DCpts == 1));
            if Changenum >= numel(Changeops)
                Changenum = numel(Changeops);
            end
            turnidx = datasample(Changeops,Changenum);
            DCpts(turnidx) = 1;
            
            DC2much = [(DCY(i_fr,1) - sum(DCpts == 1)/pts) > 0.01,...
                       (DCY(i_fr,2) - sum(DCpts == 2)/pts) > 0.01,...
                       (DCY(i_fr,3) - sum(DCpts == 3)/pts) > 0.01];
            totalturn = [totalturn; turnidx];
        end
        
        if DC2much(2) % Too little with DC2
            if      ~DC2much(1) && DC2much(3)
                Changeops = find(DCpts == 1);
            elseif  DC2much(1) && ~DC2much(3)
                Changeops = find(DCpts == 3);
            elseif  ~DC2much(1) && ~DC2much(3)
                Changeops = find(DCpts == 1 | DCpts == 3);
            end
            Changenum = round(DCY(i_fr,2)*pts - sum(DCpts == 2));
            if Changenum >= numel(Changeops)
                Changenum = numel(Changeops);
            end
            turnidx = datasample(Changeops,Changenum);
            DCpts(turnidx) = 2;
            
            DC2much = [(DCY(i_fr,1) - sum(DCpts == 1)/pts) > 0.01,...
                       (DCY(i_fr,2) - sum(DCpts == 2)/pts) > 0.01,...
                       (DCY(i_fr,3) - sum(DCpts == 3)/pts) > 0.01];
            totalturn = [totalturn; turnidx];
        end
        
        if DC2much(3) % Too little with DC3
            if      ~DC2much(1) && DC2much(2)
                Changeops = find(DCpts == 1);
            elseif  DC2much(1) && ~DC2much(2)
                Changeops = find(DCpts == 2);
            elseif  ~DC2much(1) && ~DC2much(2)
                Changeops = find(DCpts == 1 | DCpts == 2);
            end
            Changenum = round(DCY(i_fr,3)*pts - sum(DCpts == 3));
            if Changenum >= numel(Changeops)
                Changenum = numel(Changeops);
            end 
            turnidx = datasample(Changeops,Changenum);
            DCpts(turnidx) = 3;
            
            DC2much = [(DCY(i_fr,1) - sum(DCpts == 1)/pts) > 0.01,...
                       (DCY(i_fr,2) - sum(DCpts == 2)/pts) > 0.01,...
                       (DCY(i_fr,3) - sum(DCpts == 3)/pts) > 0.01];
            totalturn = [totalturn; turnidx];
        end
        
%         totalturn = unique(totalturn);
%         DCchange = false(pts,1);
%         DCchange(totalturn) = true;
%         
%         % Find all spots travelled long enough for a change in direction
%         % and reset counter
%         Dchangeidx = T0 > Framesec;
%         T0 = T0 + 1;
%         T0(Dchangeidx) = T0(Dchangeidx) - Framesec;
%         
%         % Find spots with change in direction or DC, with both DC
%         Changepts1 = (DCpts == 1) & Allchange;
%         Changepts2 = (DCpts == 2) & Allchange;
%         Changepts3 = (DCpts == 3) & Allchange;

        Changepts1 = (DCpts == 1);
        Changepts2 = (DCpts == 2);
        Changepts3 = (DCpts == 3);
        d1pts = sum(Changepts1);
        d2pts = sum(Changepts2);
        d3pts = sum(Changepts3);
        
        % Get new change in X and Y if neccessary
        XYC(Changepts1,:) = 2*(randi([0,1],d1pts,2)-1/2) .* normrnd(Dac(1),SigmaD(1),d1pts,2);
        XYC(Changepts2,:) = 2*(randi([0,1],d2pts,2)-1/2) .* normrnd(Dac(2),SigmaD(2),d2pts,2);
        XYC(Changepts3,:) = 2*(randi([0,1],d3pts,2)-1/2) .* normrnd(Dac(3),SigmaD(3),d3pts,2);
        
        % Apply change to positions
        X0 = X0 + XYC(:,1);
        Y0 = Y0 + XYC(:,2);
        
        % Find spots that have gone out of bounds
        Xstay = X0<0 | X0>TLx;
        Ystay = Y0<0 | Y0>Ly;
        
        ZeroStay=0;
        if isempty(Xstay) || Xstay(1)==0
            ZeroStay=ZeroStay+1;
        else
        % Get new change in X and Y
        XYC(Xstay,1) = -XYC(Xstay,1);
        XYC(Ystay,2) = -XYC(Ystay,2);
        
        % Limit Movement of out-of-bound spots
        X0(Xstay) = X0(Xstay) + 2*XYC(Xstay,1);
        Y0(Ystay) = Y0(Ystay) + 2*XYC(Ystay,2);
        end
        
        clear Xstay Ystay
        
        % Change in Intensity
        Ichange = normrnd(ini.MeanI,ini.SigmaI,pts,1);
        I0 = I0.*Ichange;

        %% Merging
        
        % XYC,Dac,SigmaD T0, Framesec
        
        if ini.tog_merging == 1;
            [X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts,npts,NewL] = ...
                Tig_merging(X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts,npts,NewL,ini,i_fr);
        end

        %% Splitting
        if ini.tog_split == 1;
            [X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts,npts,NewL] = ...
                Tig_splitting(X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts,npts,NewL,ini,i_fr);
        end
        
        %% Spot Disappearence
        if ini.tog_spotdis == 1;
            [X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts] = ...
                Tig_Dissappearance(X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts,ini);
        end

        %% Spot Appearence
        if ini.tog_spotapp == 1; 
            [X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts,npts,NewL] = ...
                Tig_Appearance(X0,Y0,Z0,XYC,T0,I0,L0,DCpts,npts,meanI,NewL,ini,TLx,Ly,Lz);
        end

        %% Remove low intensity spots      
        [X0,Y0,Z0,XYC,T0,I0,L0,DCpts,pts] = ...
            Tig_IntensityBound(X0,Y0,Z0,XYC,T0,I0,L0,DCpts,ini);
        
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
            D0(:,3) = D0(:,3) + tif;

            blinkaction = find(D0(:,3)>D0(:,4));
            ini.Nblink = numel(blinkaction);

            % Blinkaction
            for i_blink = blinkaction';
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
        % Compute real intensity (with blinking) (using histc to take into
        % account the appearence and dissapearence of spots
        J0 = I0;
        if ini.tog_blinking == 1;
            if ~isempty(D0)

                blinks = D0(:,1);
                blinkspots = zeros(numel(blinks),1);
                for ti = 1:numel(blinks)
                    blinkspots(ti) = find(L0==blinks(ti));
                end

                blinkI = I0(blinkspots).*D0(:,2);
                J0(blinkspots) = blinkI;        
            end
        end
        
        clear NewL sortd 
    end
    
    % Choose to output file or plot simulated paths
    if ini.tog_plotpath == 0;
    
        %create file name
        tmp1=clock;
        fname=['BlurLab_Tigercreate_' date '_' num2str(tmp1(4)) '-' num2str(tmp1(5)) '-' num2str(tmp1(6))...
            '_pts-' num2str(pts) '_meanI-' num2str(meanI) '_D-' num2str(DC(:)') '_Lx-' num2str(Lx) '_Ly-' ...
            num2str(Ly) '_Lz-' num2str(Lz) '.txt'];

        % write output file
        oldfolder = pwd;
        ini.Tigfolder = uigetdir(pwd);
        cd(ini.Tigfolder)
        BlurLab_text(Xout,Yout,Zout,Iout,Fout,Lout,fname)
        cd(oldfolder);
    
    else
            
        scatvar = 0.5;
        Nspots = max(Lout);
        Nentry = numel(Lout);
        Paths = cell(Nspots,4);

        for entry = 1:Nentry
            thisL = Lout(entry);
            thisx = Paths{thisL,1};
            thisy = Paths{thisL,2};
            thisi = Paths{thisL,3};
            thisf = Paths{thisL,4};

            thisx = [thisx, Xout(entry)];
            thisy = [thisy, Yout(entry)];
            thisi = [thisi, Iout(entry)];
            thisf = [thisf, Fout(entry)];
            
            thisi(thisi==0) = 0.001;

            Paths{thisL,1} = thisx;
            Paths{thisL,2} = thisy;
            Paths{thisL,3} = thisi;
            Paths{thisL,4} = thisf;

            clear thisx thisy thisi thisf  
        end

        figure
        axis([0, max(Xout), 0, max(Yout)])
        hold on
        for ploti = 1:Nspots
            thisc = rand(1,3);
            plot(Paths{ploti,1},Paths{ploti,2},'Color',thisc)
            scatter(Paths{ploti,1},Paths{ploti,2},scatvar*Paths{ploti,3},thisc)
        end
        hold off
    end
end



