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

SigmaR = 0.3;               % Variance for normal dis. used for change in R (displacement)
SigmaA = 0.3;               % Variance for normal dis. used for change in A (angle)
SigmaI = 0.2;               % Variance for normal dis. used for change in I (Intensity)
MergeD = mean(Lx,Ly)/20;    % Distance between two spots needed for merging
CMerge = 0.2;               % Chance for a merging event (if within distance)
CSplit = 0.05;              % Chance for a splitting event
DSplit = 0.5;               % Change in displacement for one particle
CDis = 0.02;                % Chhance for a disappearence event
Apois = 0.5;                % Number of spots appearences (poisson distribution)
Nblink = 9;                 % Number of blinking spots
Tblink = 3;                 % Average blinking duration (in frames)
Cblink1 = 0.05;             % Chance of becoming a blinking spot (new spots)
Cblink2 = 0.01;             % Chance of becoming a stable spot (blinking spots)

tog_merging = 1;
tog_split = 1;
tog_spotdis = 1;
tog_spotapp = 1;
tog_blinking = 1;

%%

D0(:,1) = round(pts*rand(Nblink,1));
D0(:,2) = round(rand(Nblink,1));
D0(:,3) = round(Tblink*rand(Nblink,1));
D0(:,4) = raylrnd(Tblink,Nblink,1);
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
    J0 = I0;
    J0(blinkspots) = blinkI;
    
    %% Merging
    if tog_merging == 1;
        [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
    merging(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,MergeD,CMerge,i_fr);
    end
    
    %% Splitting
    if tog_split == 1;
        [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
    splitting(X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL,CSplit,DSplit,i_fr,meanI);

    end
    
    %% Spot Disappearence
    if tog_spotdis == 1;
        [X0,Y0,Z0,I0,L0,A0,R0,pts] = ...
            Dissappearance(X0,Y0,Z0,I0,L0,A0,R0,pts,CDis);
    end
    
    %% Spot Appearence
    if tog_spotapp == 1; 
        [X0,Y0,Z0,I0,L0,A0,R0,pts,npts,NewL] = ...
            Appearance(X0,Y0,Z0,I0,L0,A0,R0,npts,meanI,NewL,Apois,Lx,Ly,Lz,D);
    end
    
    %% Blinking
    if tog_blinking == 1;
        
        % Remove gone spots
        D0(find(ismember(D0(:,1),L0)==0),:) = [];
        
        % Add new blinking spots
        Nnspots = length(NewL);
        if ~Nnspots == 0;
            Newidx = find(rand(Nnspots,1) < Cblink1);
            Dnew(:,1) = NewL(Newidx);
            Dnew(:,2) = round(rand(numel(Newidx),1));
            Dnew(:,3) = round(Tblink*rand(numel(Newidx),1));
            Dnew(:,4) = raylrnd(Tblink,numel(Newidx),1);
            
            D0 = [D0; Dnew];
        end
        
        % Add frame to blinkcounter
        D0(:,3) = D0(:,3) + 1;
        
        blinkaction = find(D0(:,3)>D0(:,4));
        Nblink = numel(blinkaction);
        
        % Blinkaction
        for i_blink = blinkaction;
            if D0(i_blink,2) == 1;
                D0(i_blink,2) = 0;
            elseif D0(i_blink,2) == 0;
                D0(i_blink,2) = 1;
            end
            D0(i_blink,3) = 0;
            D0(i_blink,4) = raylrnd(Tblink);
        end
        
        % Remove some blinking spots
        Nspots = size(D0,1);
        Ridx = find(rand(Nspots,1) < Cblink2);
        D0(Ridx,:) = [];
        
        % sort blinking spots        
        [~, sortd] = sort(D0(:,1));
        D0 = D0(sortd,:);
    end

    
    Rchange = normrnd(1,SigmaR,pts,1);
    Achange = normrnd(1,SigmaA,pts,1);
    Ichange = normrnd(1,SigmaI,pts,1);
    
    R0 = R0.*Rchange;
    A0 = A0.*Achange;
    I0 = I0.*Ichange;
    

    clear NewL sortd 
end
    



    %I know this isn't quite right
%     X0(X0<0)=0;
%     X0(X0>Lx)=Lx;
%     Y0(Y0<0)=0;
%     Y0(Y0>Ly)=Ly;
%     Z0(Z0<0)=0;
%     Z0(Z0>Lz)=Lz;

%write output file
% BlurLab_text(Xout,Yout,Zout,Iout,Fout,Lout,fname)




    





