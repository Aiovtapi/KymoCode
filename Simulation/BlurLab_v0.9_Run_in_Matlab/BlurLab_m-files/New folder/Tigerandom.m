%Mark Shui Hu
%Random 3D Diffusion of particles
clc
clear all

nframes = 5;
pts = 50;
meanI = 20;
D = 1;
Lx = 50;
Ly = 50;
Lz = 1;


% function random_box2(nframes,pts,meanI,D,Lx,Ly,Lz)
%create file name
tmp1=clock;
fname=['BlurLab_rand_output_' date '_' num2str(tmp1(4)) '-' num2str(tmp1(5)) '-' num2str(tmp1(6))...
    '_pts-' num2str(pts) '_meanI-' num2str(meanI) '_D-' num2str(D) '_Lx-' num2str(Lx) '_Ly-' num2str(Ly) '_Lz-' num2str(Lz) '.txt'];

ini.SigmaR = 0.3;               % Variance for normal dis. used for change in R (displacement)
ini.SigmaA = 0.3;               % Variance for normal dis. used for change in A (angle)
ini.SigmaI = 0.2;               % Variance for normal dis. used for change in I (Intensity)
ini.MergeD = mean(Lx,Ly)/20;    % Distance between two spots needed for merging
ini.CMerge = 0.2;               % Chance for a merging event (if within distance)
ini.CSplit = 0.05;              % Chance for a splitting event
ini.DSplit = 0.5;               % Change in displacement for one particle
ini.CDis = 0.02;                % Chhance for a disappearence event
ini.Apois = 0.5;                % Number of spots appearences (poisson distribution)
ini.Nblink = 9;                 % Number of blinking spots
ini.Tblink = 3;                 % Average blinking duration (in frames)
ini.Cblink1 = 0.05;             % Chance of becoming a blinking spot (new spots)
ini.Cblink2 = 0.01;             % Chance of becoming a stable spot (blinking spots)

ini.tog_merging = 1;
ini.tog_split = 1;
ini.tog_spotdis = 1;
ini.tog_spotapp = 1;
ini.tog_blinking = 1;

%%

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
