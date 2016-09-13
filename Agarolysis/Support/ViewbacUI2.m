function [skip,fault,previous,Rspot] = ViewbacUI2(init,chans,f,Bacpics,Bacmesh,X,BX,celli,Title)

frames = size(Bacmesh{1},2);
[~, cells] = size(X);
skip = 0;
fault = 0;
previous = 0;
Rspot = [];
Nchans = numel(chans);

currentbac = ['Bacpic ',num2str(celli),'/',num2str(cells)];

%% Create push buttons and indicators
Skip = uicontrol('Style', 'pushbutton', 'String', 'Skip all',...
    'Position', [10 10 50 30],...
    'Callback',@Skipall);    

bac = uicontrol('Style','text',...
    'Position',[150 15 120 15],...
    'String',currentbac);

Prev = uicontrol('Style', 'pushbutton', 'String', 'Previous Cell',...
    'Position', [280 10 100 30],...
    'Callback', @GoPrev);

Disr = uicontrol('Style', 'pushbutton', 'String', 'Disregard Cell',...
    'Position', [400 10 100 30],'BackgroundColor',[1 0.6 0.6],...
    'Callback', @Savefault);

Next = uicontrol('Style', 'pushbutton', 'String', 'Next Cell',...
    'Position', [520 10 100 30],...
    'Callback',@Closefigure);

clcr = uicontrol('Style','text',...
    'Position',[1000 15 100 15],...
    'String','Click spot to remove');

clr = uicontrol('Style', 'pushbutton', 'String', 'Clear clicks',...
    'Position', [1110 10 80 30],...
    'Callback',@Clearclick);    

%% Plot first frame of bacpic, mesh and spots
for chan = chans;
    ax(Nchans) = subplot(1,Nchans,chan);
    set(gca,'tag',num2str(chan))

    title(Title{chan})

    x = X{chan,celli};
    bx = BX{chan,celli};

    hold on
    imagesc(Bacpics{chan}{celli,1});
    plot(Bacmesh{chan}{celli,1}(:,1),Bacmesh{chan}{celli,1}(:,2),'w',...
        Bacmesh{chan}{celli,1}(:,3),Bacmesh{chan}{celli,1}(:,4),'w','LineWidth',2)
    for spoti = 1:length(x)
        plot(x{spoti}(1,2),x{spoti}(1,4),'rx','LineWidth',2)
    end
    for spoti = 1:length(bx)
        plot(bx{spoti}(1,2),bx{spoti}(1,4),'kx','LineWidth',2)
    end
    axis off
    hold off
    clear x bx
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

set(f,'HitTest','off')                          % Necessary for clicking on plot
set(f,'WindowButtonDownFcn',@clicky)            % Action for clicks

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
            title(Title{chan})
            
            x = X{chan,celli};
            bx = BX{chan,celli};
            
            hold on                             % Plot the bacpic, mesh and spots
            imagesc(Bacpics{chan}{celli,frami})
            plot(Bacmesh{chan}{celli,frami}(:,1),...
                Bacmesh{chan}{celli,frami}(:,2),'w',...
                Bacmesh{chan}{celli,frami}(:,3),...
                Bacmesh{chan}{celli,frami}(:,4),'w','LineWidth',2)
            for spoti = 1:length(x)
                plot(x{spoti}(frami,2),x{spoti}(frami,4),'rx','LineWidth',2)
            end
            for spoti = 1:length(bx)
                plot(bx{spoti}(frami,2),bx{spoti}(frami,4),'kx','LineWidth',2)
            end
            axis off
            hold off
            clear x bx
        end
        
        for respot = 1:size(Rspot,1)                    % plot clicked spots
            thisx = BX{Rspot(respot,1),celli}{Rspot(respot,2)}(frami,:);
            subplot(1,Nchans,Rspot(respot,1))
            hold on
            plot(thisx(2),thisx(4),'rx','LineWidth',2)
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

    function clicky(gcbo,eventdata,handles)
        clickXY = get(gca,'CurrentPoint');              % Get click values
        clickxy = [clickXY(1,1),clickXY(1,2)];          % xy of clicked point
        clickchan = str2double((get(gca,'tag')));       % channel clicked
        
        if frames > 1
            frami = round(sld.Value);                   % get current frame
        else
            frami = 1;
        end
        
        thisx = BX{clickchan,celli};                     % get spots information
        spots = size(thisx,2);
        
        spotx = zeros(spots,2);
        for spotn = 1:spots
            spotx(spotn,1) = thisx{spotn}(frami,2);     % prepare spots coordinates
            spotx(spotn,2) = thisx{spotn}(frami,4);
        end
        
        
        [minval,minidx] = ...                           % Find minimal spot distance
            min(sqrt(sum(bsxfun(@minus,clickxy,spotx).^2,2)));
        
        % If click close position to a spot, the spot will be registered.
        % If not, a 
        if minval < 0.5
            removespot = [clickchan,minidx];
            hold on
            plot(spotx(minidx,1),spotx(minidx,2),'rx','LineWidth',2);
            hold off
        else
            removespot = [];
            bound = 2; 
            Rxy = round(clickxy);
            
            if ~any(Rxy<1)
                thisbac = Bacpics{chan}{celli,frami};
                padbac = padarray(thisbac,[bound,bound]);
                baccrop = [Rxy, 2*bound, 2*bound];
                croppedbac = imcrop(padbac, baccrop);
%                 figure(2)
%                 imagesc(croppedbac)
                x = GaussFitSimedit_ViewBacUI2(init,chan,croppedbac,init.IPTP(chan))
            end
        end
        
        Rspot = [Rspot; removespot];
        Rspot = unique(Rspot, 'rows');
    end

    function Clearclick(hObject, eventdata, handles)
        
        if frames > 1
            frame = round(sld.Value);                   % get current frame
        else
            frame = 1;
        end
        
        for respot = 1:size(Rspot,1);                   % replot clicked spot
            thisx = BX{Rspot(respot,1),celli}{Rspot(respot,2)}(frame,:);
            subplot(1,Nchans,Rspot(respot,1))
            hold on
            plot(thisx(2),thisx(4),'kx','LineWidth',2)
            hold off
        end
        
        Rspot = [];
    end
end


