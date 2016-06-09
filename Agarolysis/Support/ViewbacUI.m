function [skip,fault] = ViewbacUI(Bacpics,Mesh,x,bx,celli)

    f = figure('Visible','on');
    title(['Review bacpic no. ',num2str(celli)])
    
    frames = size(Mesh,2);
    skip = 0;
    fault = 0;
       
%% Create push buttons

    Disr = uicontrol('Style', 'pushbutton', 'String', 'Disregard Cell',...
        'Position', [100 10 100 30],...
        'Callback', @Savefault);
    
    Next = uicontrol('Style', 'pushbutton', 'String', 'Next Cell',...
        'Position', [220 10 100 30],...
        'Callback',@Closefigure);
    
    Skip = uicontrol('Style', 'pushbutton', 'String', 'Skip all',...
        'Position', [10 10 50 30],...
        'Callback',@Skipall);    

%% Plot first frame

    hold on
    imagesc(Bacpics{celli,1});
    plot(Mesh{celli,1}(:,1),Mesh{celli,1}(:,2),'w',...
        Mesh{celli,1}(:,3),Mesh{celli,1}(:,4),'w','LineWidth',2)
    for spoti = 1:length(x)
        plot(x{spoti}(1,2),x{spoti}(1,4),'rx','LineWidth',2)
    end
    for spoti = 1:length(bx)
        plot(bx{spoti}(1,2),bx{spoti}(1,4),'kx','LineWidth',2)
    end
    axis off
    hold off
    
%% Add slider if theres more than 1 frame        
	if frames > 1
        sld = uicontrol('Style', 'slider',...
            'Min',1,'Max',frames,'Value',1,...
            'Position', [340 12 200 18],...
            'Callback', @selectframe); 
    
        txt = uicontrol('Style','text',...
            'Position',[380 30 120 15],...
            'String','Frame');
	end
    
    uiwait(f) % Wait for button click

%% functions for the buttons and sliders

    function selectframe(source,callbackdata)
        frami = round(source.Value);
        
        hold on
        imagesc(Bacpics{celli,frami})
        plot(Mesh{celli,frami}(:,1),Mesh{celli,frami}(:,2),'w',...
            Mesh{celli,frami}(:,3),Mesh{celli,frami}(:,4),'w','LineWidth',2)
        for spoti = 1:length(x)
            plot(x{spoti}(frami,2),x{spoti}(frami,4),'rx','LineWidth',2)
        end
        for spoti = 1:length(bx)
            plot(bx{spoti}(frami,2),bx{spoti}(frami,4),'kx','LineWidth',2)
        end
        axis off
        hold off
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
end