function memupdate_Callback(hObject, eventdata, handles)
%get and set memory
try
    [meminfo,~]=memory;
    set(handles.text_totmem,'String',num2str(1/10*round(meminfo.MemAvailableAllArrays/1e5)));
    set(handles.text_availmem,'String',num2str(1/10*round(meminfo.MaxPossibleArrayBytes/1e5)));
catch err1
    set(handles.text_totmem,'String','na');
    set(handles.text_availmem,'String','na');
end