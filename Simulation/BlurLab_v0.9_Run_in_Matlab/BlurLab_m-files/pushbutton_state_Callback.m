function pushbutton_state_Callback(hObject, eventdata, handles)

guidata(hObject,handles)

ans1=questdlg('Are you sure you want to save the current program state as default?', ...
    'Set Defaults','Yes', 'No', 'Not as Default', 'Yes');

if or(strcmp('Yes',ans1),strcmp('Not as Default',ans1))
    
    if strcmp('Yes',ans1)
        l1=dir('BlurLab_Prefs.txt');
        if ~isempty(l1)
            delete BlurLab_Prefs.txt
        end
        fid=fopen('BlurLab_Prefs.txt','w+');
    else
        [fname_def,cd_def]=uiputfile(['BlurLab_settings_' date '.txt'],'Save a setting preference file ...');
        
        if fname_def==0
            return
        else
            fid=fopen([cd_def fname_def],'w+');
        end
    end
    
    npars=60;
    pars=zeros(1,npars);
    
    %pane 1
    pars(1)=str2num(get(handles.edit_startim,'String')); % start frame
    pars(2)=str2num(get(handles.edit_endim,'String')); % end frame
    pars(3)=get(handles.checkbox_boxcar,'Value'); % performs boxcar
    pars(4)=str2num(get(handles.edit_boxcar,'String')); % boxcar size
    pars(5)=get(handles.checkbox_histogram,'Value'); % show histogram
    pars(6)=get(handles.checkbox_log,'Value'); % log
    pars(7)=get(handles.checkbox_linear,'Value'); % linear
    pars(8)=get(handles.checkbox_disk,'Value'); % memory uses disk
    pars(9)=get(handles.checkbox_ram,'Value'); % memory uses ram
    
    %pane 2
    pars(10)=str2num(get(handles.edit_impxlsize,'String')); % image calibration
    pars(11)=str2num(get(handles.edit_psfdz,'String')); % psf dz
    pars(12)=str2num(get(handles.edit_outY,'String')); % output size Y
    pars(13)=str2num(get(handles.edit_outX,'String')); % output size X
    pars(14)=str2num(get(handles.edit_offsetx,'String')); % offget X
    pars(15)=str2num(get(handles.edit_offsety,'String')); % offget Y
    pars(16)=get(handles.checkbox_flip,'Value'); % flip image
    pars(17)=str2num(get(handles.edit_border,'String')); %image border
    
    %pane 3
    pars(18)=str2num(get(handles.edit_sfc,'String')); % scale factor
    pars(19)=get(handles.checkbox_psfnorm,'Value'); % normalize psf
    pars(20)=get(handles.checkbox_noise,'Value'); % gauss noise
    pars(21)=str2num(get(handles.edit_mean,'String')); % mean
    pars(22)=str2num(get(handles.edit_std,'String')); % std
    pars(23)=get(handles.checkbox_poiss,'Value'); % poiss noise
    pars(24)=str2num(get(handles.edit_poiss,'String')); % multiplier
    pars(25)=str2num(get(handles.edit_basal,'String')); % basal poisson
    pars(26)=get(handles.checkbox_bleach,'Value'); % photo-bleach
    pars(27)=get(handles.checkbox_stochastic,'Value'); % stochastic bleaching
    pars(28)=str2num(get(handles.edit_tau,'String')); % tau
    pars(29)=get(handles.checkbox_usenoise,'Value'); % playback
    
    %pane 4
    pars(30)=str2num(get(handles.edit_zslice,'String')); % z slice
    pars(31)=get(handles.checkbox_simZ,'Value'); % simulate z stack
    pars(32)=str2num(get(handles.edit_minz,'String')); % zmin
    pars(33)=str2num(get(handles.edit_maxz,'String')); % zmax
    pars(34)=str2num(get(handles.edit_dzstack,'String')); % dz
    
    %pane 5
    pars(35)=get(handles.checkbox_write,'Value'); % save images
    pars(36)=get(handles.checkbox_savepars,'Value'); % save parameters with images
    pars(37)=get(handles.checkbox_8bit,'Value'); % 8-bit
    pars(38)=get(handles.checkbox_16bit,'Value'); % 16-bit
    
    %pane 6
    pars(39)=get(handles.checkbox_frap,'Value'); % FRAP setting
    pars(40)=str2num(get(handles.edit_frapframe,'String')); % FRAP frame start
    pars(41)=str2num(get(handles.edit_frap,'String')); % FRAP percentage
    pars(42)=get(handles.checkbox_tirf,'Value'); % TIRF setting
    pars(43)=str2num(get(handles.edit_tirfdecay,'String')); % TIRF decay length
    pars(44)=str2num(get(handles.edit_tirfZ,'String')); % TIR Z plane location
    
    %pane  7
    pars(45)=str2num(get(handles.edit_points,'String')); % points
    pars(46)=str2num(get(handles.edit_meanint,'String')); % mean intensity
    pars(47)=str2num(get(handles.edit_movestd,'String')); % diffusion coeff.
    pars(48)=str2num(get(handles.edit_randframes,'String')); % frames
    pars(49)=str2num(get(handles.edit_Lx,'String')); % X box size
    pars(50)=str2num(get(handles.edit_Ly,'String')); % Y box size
    pars(51)=str2num(get(handles.edit_Lz,'String')); % Z box size
    
    %figure pane
    pars(52)=get(handles.checkbox_hot,'Value'); % hot
    pars(53)=get(handles.checkbox_jet,'Value'); % jet
    pars(54)=get(handles.checkbox_gray,'Value'); % gray
    pars(55)=get(handles.checkbox_fix,'Value'); % fix contrast
    pars(56)=str2num(get(handles.edit_lowfix,'String')); % int min
    pars(57)=str2num(get(handles.edit_highfix,'String')); % int max
    
    %playback pane
    pars(58)=get(handles.checkbox_loop,'Value'); %loop playback
    pars(59)=str2num(get(handles.edit_fps,'String')); % frames per second
    
    %misc
    pars(60)=get(handles.popupmenu1,'Value'); % menu choice
    
    %setup output preferences file
    for i=1:npars
        fprintf(fid,'%6.4f\n',pars(i));
    end
    fclose(fid);
end




