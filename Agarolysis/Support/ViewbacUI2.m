function [skip,fault,previous] = ViewbacUI2(Bacpics,Mesh,X,BX,celli,Title)

    frames = size(Mesh{1},2);
    [chans, cells] = size(X);
    skip = 0;
    fault = 0;
    previous = 0;
    
    f = figure('Visible','on','Position',[350 350 1200, 400]);
    currentbac = ['Bacpic ',num2str(celli),'/',num2str(cells)];
       
%% Create push buttons
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

%% Plot first frame
    for chan = 1:chans;
        subplot(1,chans,chan)
        title(Title{chan})
        
        x = X{chan,celli};
        bx = BX{chan,celli};
        
        hold on
        imagesc(Bacpics{chan}{celli,1});
        plot(Mesh{chan}{celli,1}(:,1),Mesh{chan}{celli,1}(:,2),'w',...
            Mesh{chan}{celli,1}(:,3),Mesh{chan}{celli,1}(:,4),'w','LineWidth',2)
        for spoti = 1:length(x)
            plot(x{spoti}(1,2),x{spoti}(1,4),'rx','LineWidth',2)
        end
        for spoti = 1:length(bx)
            plot(bx{spoti}(1,2),bx{spoti}(1,4),'kx','LineWidth',2)
        end
        axis off
        hold off
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
    
    uiwait(f) % Wait for button click

%% functions for the buttons and sliders

    function selectframe(source,callbackdata)
        frami = round(source.Value);
        currentframe = ['Frame ',num2str(frami),'/',num2str(frames)];
        
        frm = uicontrol('Style','text',...
            'Position',[900 15 50 15],...
            'String',currentframe);

        for chan = 1:chans;
            subplot(1,chans,chan)
            title(Title{chan})
            
            x = X{chan,celli};
            bx = BX{chan,celli};
            
            hold on
            imagesc(Bacpics{chan}{celli,frami})
            plot(Mesh{chan}{celli,frami}(:,1),Mesh{chan}{celli,frami}(:,2),'w',...
                Mesh{chan}{celli,frami}(:,3),Mesh{chan}{celli,frami}(:,4),'w','LineWidth',2)
            for spoti = 1:length(X{chan,celli})
                plot(x{spoti}(frami,2),x{spoti}(frami,4),'rx','LineWidth',2)
            end
            for spoti = 1:length(bx)
                plot(bx{spoti}(frami,2),bx{spoti}(frami,4),'kx','LineWidth',2)
            end
            axis off
            hold off
        end
    end

    function Closefigure(hObject, eventdata, handles)
        close(f)
    end

    function Skipall(hObject, eventdata, handles)
        skip = 1;
        close(f)
    end

    function Savefault(hObject, eventdata, handles)
        fault = 1;
        close(f);
    end

    function GoPrev(hObject, eventdata, handles)
        previous = 1;
        close(f);
    end
end