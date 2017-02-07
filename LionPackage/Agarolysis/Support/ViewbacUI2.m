function [skip,fault,previous,Rspot,Nspot] = ViewbacUI2(init,chans,f,Bacpics,Bacmesh,CBacmask,DataStruct,cells,celli)

frames = size(Bacmesh{1},2);
skip = 0;
fault = 0;
previous = 0;
Rspot = [];
Nspot = [];
Nchans = numel(chans);

Xsize = 10;
Lwidth = 3;

currentbac = ['Bacpic ',num2str(celli),'/',num2str(cells)];

%% Create push buttons and indicators
Skip = uicontrol('Style', 'pushbutton', 'String', 'Skip all (S)',...
    'Position', [10 10 100 30],...
    'Callback',@Skipall);    

bac = uicontrol('Style','text',...
    'Position',[150 15 120 15],...
    'String',currentbac);

Prev = uicontrol('Style', 'pushbutton', 'String', 'Previous Cell (J)',...
    'Position', [280 10 100 30],...
    'Callback', @GoPrev);

Disr = uicontrol('Style', 'pushbutton', 'String', 'Disregard Cell (D)',...
    'Position', [400 10 100 30],'BackgroundColor',[1 0.6 0.6],...
    'Callback', @Savefault);

Next = uicontrol('Style', 'pushbutton', 'String', 'Next Cell (K)',...
    'Position', [520 10 100 30],...
    'Callback',@Closefigure);

clcr = uicontrol('Style','text',...
    'Position',[1000 15 100 15],...
    'String','Click spot to remove');

clr = uicontrol('Style', 'pushbutton', 'String', 'Clear clicks (C)',...
    'Position', [1110 10 80 30],...
    'Callback',@Clearclick);

%% Plot first frame of bacpic, mesh and spots
for chan = chans;
    ax(Nchans) = subplot(1,Nchans,chan);   
    
    set(gca,'tag',num2str(chan))

    title(init.channels{chan})
    
    x = DataStruct(chan,celli).x;
    bx = DataStruct(chan,celli).bx;

    hold on
    imagesc(Bacpics{chan}{celli,1});
    plot(Bacmesh{chan}{celli,1}(:,1),Bacmesh{chan}{celli,1}(:,2),'w',...
        Bacmesh{chan}{celli,1}(:,3),Bacmesh{chan}{celli,1}(:,4),'w','LineWidth',Lwidth)
%     for spoti = 1:length(x)
%         plot(x{spoti}(1,2),x{spoti}(1,4),'rx','LineWidth',2)
%     end
    for spoti = 1:length(bx)
        plot(bx{spoti}(1,2),bx{spoti}(1,4),'kx','LineWidth',Lwidth,'MarkerSize',Xsize)
    end
    axis off
    hold off
    clear x bx
    
    pbound = 10;
    fsize = get(gcf,'Position');
    FXsize = fsize(3);
    FYsize = fsize(4);
    Fratio = FXsize/Nchans/FYsize;
    
    XALim = get(gca,'XLim');
    YALim = get(gca,'YLim');
    XAsize = abs(XALim(2)-XALim(1));
    YAsize = abs(YALim(2)-YALim(1));
    
    if XAsize>YAsize*Fratio
        pxpos = FXsize/Nchans*(chan-1) + pbound;
        xlength = round((FXsize-(Nchans*pbound))/Nchans);
        ylength = round(YAsize*Fratio/XAsize*xlength);
        pypos = round((FYsize-ylength+50)/2);       
    else
        ylength = round(0.75*FYsize);
        pypos = 50;
        xlength = round(XAsize/(YAsize*Fratio)*ylength);
        pxpos = FXsize/Nchans*(chan-1) + round((FXsize/Nchans-xlength)/2);
    end
    POS = [pxpos/FXsize, pypos/FYsize, xlength/FXsize, ylength/FYsize];
    set(gca,'Position',POS)
end

%% Add slider if theres more than 1 frame        
if frames > 1
    sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',frames,'Value',1,...
        'Position', [640 12 250 18],...
        'Callback', @selectframe); 

    txt = uicontrol('Style','text',...
        'Position',[680 30 120 15],...
        'String','Frame');
end

%%
Rspots = [];

figure(f)
set(f,'HitTest','off')                          % Necessary for clicking on plot
set(f,'WindowButtonDownFcn',@clicky)            % Action for clicks
set(f,'KeyPressFcn',@shortcuts);

uiwait(f)                                       % Wait for button click


%% functions for the buttons and sliders

    function selectframe(source,callbackdata)   % Plotting for change of frame
        frami = round(source.Value);
        currentframe = ['Frame ',num2str(frami),'/',num2str(frames)];
        
    frm = uicontrol('Style','text',...          % Add frame counter
            'Position',[900 15 50 15],...
            'String',currentframe);

        for chan = chans;
            ax(Nchans)=subplot(1,Nchans,chan);
            title(init.channels{chan})
            
            x = DataStruct(chan,celli).x;
            bx = DataStruct(chan,celli).bx;
            
            hold on                             % Plot the bacpic, mesh and spots
            imagesc(Bacpics{chan}{celli,frami})
            plot(Bacmesh{chan}{celli,frami}(:,1),...
                Bacmesh{chan}{celli,frami}(:,2),'w',...
                Bacmesh{chan}{celli,frami}(:,3),...
                Bacmesh{chan}{celli,frami}(:,4),'w','LineWidth',Lwidth)
%             for spoti = 1:length(x)
%                 plot(x{spoti}(frami,2),x{spoti}(frami,4),'rx','LineWidth',2)
%             end
            for spoti = 1:length(bx)
                plot(bx{spoti}(frami,2),bx{spoti}(frami,4),'kx','LineWidth',Lwidth,'MarkerSize',Xsize)
            end
            axis off
            hold off
            clear x bx
        end
        
        for respot = 1:size(Rspot,1)                    % plot clicked spots

            thisx = DataStruct(Rspot(respot,1),celli).bx{Rspot(respot,2)}(frami,:);
            subplot(1,Nchans,Rspot(respot,1))
            hold on
            plot(thisx(2),thisx(4),'rx','LineWidth',Lwidth,'MarkerSize',Xsize)
            hold off
        end
            
    end

    function Closefigure(hObject, eventdata, handles)   % Close and continue
        uiresume(f)
        clf(f)
    end

    function Skipall(hObject, eventdata, handles)       % Skip all cells
        skip = 1;
        uiresume(f)
        clf(f)
    end

    function Savefault(hObject, eventdata, handles)     % Save faulty cell
        fault = 1;
        uiresume(f)
        clf(f)
    end

    function GoPrev(hObject, eventdata, handles)        % Go to previous cell
        previous = 1;
        uiresume(f);
        clf(f)
    end

    function shortcuts(gcbo,eventdata,handles)
        switch eventdata.Key
            case 'j'
                GoPrev
            case 'k'
                Closefigure
            case 'd'
                Savefault
            case 's'
                Skipall
            case 'c'
                Clearclick
        end
    end

    function clicky(gcbo,eventdata,handles)
        clickXY = get(gca,'CurrentPoint');              % Get click values
        clickxy = [clickXY(1,1),clickXY(1,2)];          % xy of clicked point
        clickchan = str2double((get(gca,'tag')));       % channel clicked
        
        if frames > 1
            frami = round(sld.Value);                   % get current frame
        else
            frami = 1;
        end

        thisx = DataStruct(clickchan,celli).bx;         % get spots information
        rspots = size(thisx,2);
        nspots = size(Nspot,1);
                
        spotx = zeros(rspots+nspots,2);
        for spotn = 1:rspots
            spotx(spotn,1) = thisx{spotn}(frami,2);     % prepare spots coordinates
            spotx(spotn,2) = thisx{spotn}(frami,4);
        end
        
        if ~nspots ==0 
            for spotn = 1:nspots
                spotx(rspots + spotn,1) = Nspot(spotn,2);
                spotx(rspots + spotn,2) = Nspot(spotn,4);
            end
        end
        
        [minval,minidx] = ...                           % Find minimal spot distance
            min(sqrt(sum(bsxfun(@minus,clickxy,spotx).^2,2)));
        
        % If click close position to a spot, the spot will be registered.
        % If not, a 
        if minval < 0.5
            
            if ~isempty(Rspot)
                thisRspots = find(Rspot(:,1)==clickchan);
                if numel(thisRspots) > 0
                    for i = 1:numel(thisRspots)
                        thisrpos(i,1) = thisx{Rspot(thisRspots(i),2)}(frami,2);
                        thisrpos(i,2) = thisx{Rspot(thisRspots(i),2)}(frami,4);
                    end
                    [oldspot,rspotidx] = ismember(thisrpos,spotx(minidx,:),'rows');
                else
                    oldspot = 0;
                end
            else
                oldspot = 0;
            end
            
            hold on
            if oldspot
                plot(spotx(minidx,1),spotx(minidx,2),'kx','LineWidth',Lwidth,'MarkerSize',Xsize);
            else
                plot(spotx(minidx,1),spotx(minidx,2),'rx','LineWidth',Lwidth,'MarkerSize',Xsize);
            end
            hold off
            
            
            if minidx <= rspots
                if oldspot
                    Rspot(thisRspots(rspotidx),:) = [];
                    removespot = [];
                else
                    removespot = [clickchan,minidx];
                end
            else
                Nspot(minidx - rspots,:) = [];
                removespot = [];
            end
            newx = [];
            
        else
            removespot = [];
            bound = 2; 
            Rxy = round(clickxy);
            
            if ~any(Rxy<1)
                thisbac = Bacpics{clickchan}{celli,frami};
                padbac = padarray(thisbac,[bound,bound]);
                baccrop = [Rxy, 2*bound, 2*bound];
                croppedbac = imcrop(padbac, baccrop);
                newx = GaussFitSimedit_ViewBacUI2(init,chan,croppedbac);
                newx = newx{1,1};
                newx(2) = newx(2) + Rxy(1) - bound - 1;
                newx(4) = newx(4) + Rxy(2) - bound - 1;
                if ~isempty(thisx)
                    newx(7) = thisx{1}(frami,7);
                else
                    newx(7) = sum(thisbac(logical(CBacmask{celli,frami})));
                end
                hold on
                plot(newx(2),newx(4),'kx','LineWidth',Lwidth,'MarkerSize',Xsize)
                hold off
                newx(9) = clickchan;
                newx(10) = frami;
            else
                newx = [];
            end
        end
        Nspot = [Nspot; newx];
        Nspot = unique(Nspot,'rows');
        
        Rspot = [Rspot; removespot];
        Rspot = unique(Rspot,'rows');
        
        clear thisrpos
    end

    function Clearclick(hObject, eventdata, handles)
        
        if frames > 1
            frami = round(sld.Value);                   % get current frame
        else
            frami = 1;
        end
        
%         for respot = 1:size(Rspot,1);                   % replot clicked spot
%             thisx = BX{Rspot(respot,1),celli}{Rspot(respot,2)}(frame,:);
%             subplot(1,Nchans,Rspot(respot,1))
%             hold on
%             plot(thisx(2),thisx(4),'kx','LineWidth',2)
%             hold off
%         end
        
        for chan = chans;
            ax(Nchans)=subplot(1,Nchans,chan);
            title(init.channels{chan})
            
            x = DataStruct(chan,celli).x;
            bx = DataStruct(chan,celli).bx;
            
            hold on                             % Plot the bacpic, mesh and spots
            imagesc(Bacpics{chan}{celli,frami})
            plot(Bacmesh{chan}{celli,frami}(:,1),...
                Bacmesh{chan}{celli,frami}(:,2),'w',...
                Bacmesh{chan}{celli,frami}(:,3),...
                Bacmesh{chan}{celli,frami}(:,4),'w','LineWidth',Lwidth)
            for spoti = 1:length(x)
                plot(x{spoti}(frami,2),x{spoti}(frami,4),'rx','LineWidth',Lwidth,'MarkerSize',Xsize)
            end
            for spoti = 1:length(bx)
                plot(bx{spoti}(frami,2),bx{spoti}(frami,4),'kx','LineWidth',Lwidth,'MarkerSize',Xsize)
            end
            axis off
            hold off
            clear x bx
        end
        
        Rspot = [];
        Nspot = [];
    end
end


