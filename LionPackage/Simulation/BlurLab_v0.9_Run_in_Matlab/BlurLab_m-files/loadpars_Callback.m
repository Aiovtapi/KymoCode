function loadpars_Callback(hObject, eventdata, handles)
try
    [parsname,cd_pars,~]=uigetfile({'*.txt','Text Input (*.txt)'},'Pick a parameter input file.');
    
    if parsname==0
        return
    end
    
    pars=load([cd_pars parsname]);
    if size(pars,1)==60
%pane 1
        set(handles.edit_startim,'String',num2str(pars(1))); % start frame
        set(handles.edit_endim,'String',num2str(pars(2))); % end frame
        set(handles.checkbox_boxcar,'Value',pars(3)); % performs boxcar
        set(handles.edit_boxcar,'String',num2str(pars(4))); % boxcar size
        set(handles.checkbox_histogram,'Value',pars(5)); % show histogram
        set(handles.checkbox_log,'Value',pars(6)); % log
        set(handles.checkbox_linear,'Value',pars(7)); % linear
        set(handles.checkbox_disk,'Value',pars(8)); % memory uses disk
        set(handles.checkbox_ram,'Value',pars(9)); % memory uses ram
        
        %pane 2
        set(handles.edit_impxlsize,'String',num2str(pars(10))); % image calibration
        set(handles.edit_psfdz,'String',num2str(pars(11))); % psf dz
        set(handles.edit_outY,'String',num2str(pars(12))); % output size Y
        set(handles.edit_outX,'String',num2str(pars(13))); % output size X
        set(handles.edit_offsetx,'String',num2str(pars(14))); % offget X
        set(handles.edit_offsety,'String',num2str(pars(15))); % offget Y
        set(handles.checkbox_flip,'Value',pars(16)); % flip image
        set(handles.edit_border,'String',num2str(pars(17))); %image border
        
        %pane 3
        set(handles.edit_sfc,'String',num2str(pars(18))); % scale factor
        set(handles.checkbox_psfnorm,'Value',pars(19)); % normalize psf
        set(handles.checkbox_noise,'Value',pars(20)); % gauss noise
        set(handles.edit_mean,'String',num2str(pars(21))); % mean
        set(handles.edit_std,'String',num2str(pars(22))); % std
        set(handles.checkbox_poiss,'Value',pars(23)); % poiss noise
        set(handles.edit_poiss,'String',num2str(pars(24))); % multiplier
        set(handles.edit_basal,'String',num2str(pars(25))); % basal poisson
        set(handles.checkbox_bleach,'Value',pars(26)); % photo-bleach
        set(handles.checkbox_stochastic,'Value',pars(27)); % stochastic bleaching
        set(handles.edit_tau,'String',num2str(pars(28))); % tau
        set(handles.checkbox_usenoise,'Value',pars(29)); % playback
        
        %pane 4
        set(handles.edit_zslice,'String',num2str(pars(30))); % z slice
        set(handles.checkbox_simZ,'Value',pars(31)); % simulate z stack
        set(handles.edit_minz,'String',num2str(pars(32))); % zmin
        set(handles.edit_maxz,'String',num2str(pars(33))); % zmax
        set(handles.edit_dzstack,'String',num2str(pars(34))); % dz
        
        %pane 5
        set(handles.checkbox_write,'Value',pars(35)); % save images
        set(handles.checkbox_savepars,'Value',pars(36)); % save parameters with images
        set(handles.checkbox_8bit,'Value',pars(37)); % 8-bit
        set(handles.checkbox_16bit,'Value',pars(38)); % 16-bit
        
        %pane 6
        set(handles.checkbox_frap,'Value',pars(39)); % FRAP setting
        set(handles.edit_frapframe,'String',num2str(pars(40))); % FRAP frame start
        set(handles.edit_frap,'String',num2str(pars(41))); % FRAP percentage
        set(handles.checkbox_tirf,'Value',pars(42)); % TIRF setting
        set(handles.edit_tirfdecay,'String',num2str(pars(43))); % TIRF decay length
        set(handles.edit_tirfZ,'String',num2str(pars(44))); % TIR Z plane location
        
        %pane  7
        set(handles.edit_points,'String',num2str(pars(45))); % points
        set(handles.edit_meanint,'String',num2str(pars(46))); % mean intensity
        set(handles.edit_movestd,'String',num2str(pars(47))); % diffusion coeff.
        set(handles.edit_randframes,'String',num2str(pars(48))); % frames
        set(handles.edit_Lx,'String',num2str(pars(49))); % X box size
        set(handles.edit_Ly,'String',num2str(pars(50))); % Y box size
        set(handles.edit_Lz,'String',num2str(pars(51))); % Z box size
        
        %figure pane
        set(handles.checkbox_hot,'Value',pars(52)); % hot
        set(handles.checkbox_jet,'Value',pars(53)); % jet
        set(handles.checkbox_gray,'Value',pars(54)); % gray
        set(handles.checkbox_fix,'Value',pars(55)); % fix contrast
        set(handles.edit_lowfix,'String',num2str(pars(56))); % int min
        set(handles.edit_highfix,'String',num2str(pars(57))); % int max
        
        %playback pane
        set(handles.checkbox_loop,'Value',pars(58)); %loop playback
        set(handles.edit_fps,'String',num2str(pars(59))); % frames per second
        
        %misc
        set(handles.popupmenu1,'Value',pars(60)); % menu choice
        
        guidata(hObject,handles)
    else
        warndlg('Parameters file improperly formatted -- using internal defaults.','Preference File Error')
        causeerror()
    end
catch
    %pane 1
    set(handles.edit_startim,'String','1'); % start frame
    set(handles.edit_endim,'String','1'); % end frame
    set(handles.checkbox_boxcar,'Value',0); % performs boxcar
    set(handles.edit_boxcar,'String','1'); % boxcar size
    set(handles.checkbox_histogram,'Value',1); % show histogram
    set(handles.checkbox_log,'Value',0); % log
    set(handles.checkbox_linear,'Value',1); % linear
    set(handles.checkbox_disk,'Value',0); % memory uses disk
    set(handles.checkbox_ram,'Value',1); % memory uses ram
    
    %pane 2
    set(handles.edit_impxlsize,'String','80'); % image calibration
    set(handles.edit_psfdz,'String','40'); % psf dz
    set(handles.edit_outY,'String','100'); % output size Y
    set(handles.edit_outX,'String','100'); % output size X
    set(handles.edit_offsetx,'String','0'); % offget X
    set(handles.edit_offsety,'String','0'); % offget Y
    set(handles.checkbox_flip,'Value',0); % flip image
    set(handles.edit_border,'String','500'); %image border
    
    %pane 3
    set(handles.edit_sfc,'String','1'); % scale factor
    set(handles.checkbox_psfnorm,'Value',1); % normalize psf
    set(handles.checkbox_noise,'Value',0); % gauss noise
    set(handles.edit_mean,'String','300'); % mean
    set(handles.edit_std,'String','15'); % std
    set(handles.checkbox_poiss,'Value',0); % poiss noise
    set(handles.edit_poiss,'String','1'); % multiplier
    set(handles.edit_basal,'String',0'); % basal poisson
    set(handles.checkbox_bleach,'Value',0); % photo-bleach
    set(handles.checkbox_stochastic,'Value',0); % stochastic bleaching
    set(handles.edit_tau,'String','50'); % tau
    set(handles.checkbox_usenoise,'Value',0); % playback
    
    %pane 4
    set(handles.edit_zslice,'String','0'); % z slice
    set(handles.checkbox_simZ,'Value',0); % simulate z stack
    set(handles.edit_minz,'String','-500'); % zmin
    set(handles.edit_maxz,'String','500'); % zmax
    set(handles.edit_dzstack,'String','10'); % dz
    
    %pane 5
    set(handles.checkbox_write,'Value',0); % save images
    set(handles.checkbox_savepars,'Value',0); % save parameters with images
    set(handles.checkbox_8bit,'Value',0); % 8-bit
    set(handles.checkbox_16bit,'Value',1); % 16-bit
    
    %pane 6
    set(handles.checkbox_frap,'Value',0); % FRAP setting
    set(handles.edit_frapframe,'String','1'); % FRAP frame start
    set(handles.edit_frap,'String','100'); % FRAP percentage
    set(handles.checkbox_tirf,'Value',0); % TIRF setting
    set(handles.edit_tirfdecay,'String','300'); % TIRF decay length
    set(handles.edit_tirfZ,'String','0'); % TIR Z plane location
    
    %pane  7
    set(handles.edit_points,'String','100'); % points
    set(handles.edit_meanint,'String','10'); % mean intensity
    set(handles.edit_movestd,'String','1'); % diffusion coeff.
    set(handles.edit_randframes,'String','100'); % frames
    set(handles.edit_Lx,'String','20'); % X box size
    set(handles.edit_Ly,'String','20'); % Y box size
    set(handles.edit_Lz,'String','5'); % Z box size
    
    %figure pane
    set(handles.checkbox_hot,'Value',1); % hot
    set(handles.checkbox_jet,'Value',0); % jet
    set(handles.checkbox_gray,'Value',0); % gray
    set(handles.checkbox_fix,'Value',0); % fix contrast
    set(handles.edit_lowfix,'String','0'); % int min
    set(handles.edit_highfix,'String','255'); % int max
    
    %playback pane
    set(handles.checkbox_loop,'Value',0); %loop playback
    set(handles.edit_fps,'String','10'); % frames per second
    
    %misc
    set(handles.popupmenu1,'Value',1); % menu choice
    
    guidata(hObject,handles)
end

