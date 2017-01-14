clear all
close all

ini.SigmaI = 0.7;   % Variance for normal dis. used for change in I (Intensity)
ini.MergeD = 0.1;   % Distance between two spots needed for merging
ini.CMerge = 0.1;   % Chance for a merging event (if within distance)
ini.CSplit = 0.01;   % Chance for a splitting event
ini.CDis = 0.1;       % Chance for a disappearence event
ini.Apois = 0.85;     % Number of spots appearences (poisson distribution)
ini.Nblink = 0;   % Percentage of blinking spots
ini.Tblink = 5;   % Average blinking duration (in frames)
ini.Cblink1 = 0; % Chance of becoming a blinking spot (new spots)
ini.Cblink2 = 0; % Chance of becoming a stable spot (blinking spots)
ini.MinI = 20;       % Minimal value of intensity
ini.MaxI = 400;       % Maximal value of intensity

ini.tog_merging = 1;
ini.tog_split = 1;
ini.tog_spotdis = 1;
ini.tog_spotapp = 1;
ini.tog_blinking = 0;
ini.tog_dimer = 1;
ini.tog_plotpath = 1;
ini.Tplot = 1;
ini.Cellgrowth = 1;

nframes = 50;
pts = 5;
meanI = 200;
tif = 1;
Lx = 40;
Ly = 10;
Lz = 1;

ini.DC = zeros(1,3);

ini.numDC = 1;
ini.DC(1) = 0.3;
ini.DCY(:,1) = 1;

%     ini.numDC = 2;
%     ini.DC(1) = str2double(handles.TDDC_edit_DC1.String);
%     ini.DC(2) = str2double(handles.TDDC_edit_DC2.String);
%     ini.DCY(:,2) = TDDC_createvec(handles.plotx1,handles.ploty1,nframes)'/100;
%     ini.DCY(:,1) = 1 - ini.DCY(:,2);

%%
Nruns = 20;
ranval = 15;
frval = 5;
gval = 0.4;
Title = {'Simuated data 1','Simulated data 2'};

for k = 1:2
    ranx = round(rand(1,Nruns)*ranval-ranval/2);
    ranf = round(rand(1,Nruns)*frval-frval/2);
    rang = (rand(1,Nruns)*gval-gval/2);
    P = []; I = []; F = []; L = [];

    for j = 1:Nruns
        lx = Lx+ranx(j);
        frames = nframes+ranf(j);
        ini.DCY = zeros(frames,3);
        F = [];  
        [Paths, Nspots] = TigercreateVP(frames,tif,pts,meanI,lx,Ly,Lz,ini);

        for i = 1:Nspots
            P = [P, Paths{i,1}/lx];
            I = [I, Paths{i,3}];
            F = [F, Paths{i,4}];
        end

        framelength = linspace(lx/(2+rang(j)),lx,frames);
        L = [L, framelength(F)];

        clear framelength F
    end

    Nv = 20;
    Font = 'Raleway Medium';
    screensize = get( groot, 'Screensize' );
    screensize(4) = screensize(4)/2;
    AW = 1.5;
    FS = 14;
    screensize(3) = screensize(3)/3*2;
    blue = [0, 0.5, 1];

    CLbound = [14 38];

    CL = 'Cell Length';
    p = 'Normalized Cell Position';
    SI = 'Spot Intensity';
    SIN = 'Spot Intensity / Number of Spots';
    NS = 'Number of Spots';

    %%
    figure(1)
    subplot(1,2,k)
    hold on
    scatter(single(L),single(P),single(I)/Nv,'m','filled');
    hold off
    xlabel(CL); ylabel(p)
    t=title(Title{k});
    axis([10 50 0 1])
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    set(t,'Fontname',Font,'FontSize',FS)
    set(gcf, 'position', screensize)

    %%
    figure(2)
    subplot(1,2,k)
    hold on
    scatter(P,I,'x','MarkerEdgeColor','m');
    myfit=polyfit(P,I,4);
    x=0:00.1:1;
    y=polyval(myfit,x);
    plot(x,y,'k','LineWidth',3)
    ylabel(SI); xlabel(p)
    t=title(Title{k});
    hold off
    % axis([0 1 -0.1 35])
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    set(t,'Fontname',Font,'FontSize',FS)
    set(gcf, 'position', screensize)

    %%

    figure(3)
    subplot(1,2,k)
    [numbin,edges] = histcounts(P,20);
    norm = max(numbin)/350;
    X = diff(edges);
    X = cumsum(X) - X(1)/2;
    hold on
    scatter(P,I,'x','MarkerEdgeColor','m');
    plot(X,numbin/norm,'k','LineWidth',3)
    hold off
    ylabel(SIN); xlabel(p);
    t=title(Title{k});
    axis([0 1 0 420])
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    set(t,'Fontname',Font,'FontSize',FS)
    set(gcf, 'position', screensize)

    %%
    clear kNum
    figure(4)
    subplot(1,2,k)
    bins = 15;
    thisedge2{1} = linspace(min(L),max(L),bins+1);
    thisedge2{2} = (0:bins)/bins;

    % CFP
    kNum(1,:) = L;
    kNum(2,:) = P;

    Heatmap = hist3(kNum','Edges',thisedge2);
    pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
    colormap(hot) % heat map
    ylabel(p); xlabel(CL)
    t=title(Title{k});
    grid on
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    set(t,'Fontname',Font,'FontSize',FS)
    set(gcf, 'position', screensize)
end