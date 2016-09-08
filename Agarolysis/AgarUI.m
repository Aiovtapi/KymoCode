function varargout = AgarUI(varargin)
% AGARUI MATLAB code for AgarUI.fig
%      AGARUI, by itself, creates a new AGARUI or raises the existing
%      singleton*.
%
%      H = AGARUI returns the handle to a new AGARUI or the handle to
%      the existing singleton*.
%
%      AGARUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AGARUI.M with the given input arguments.
%
%      AGARUI('Property','Value',...) creates a new AGARUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AgarUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AgarUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AgarUI

% Last Modified by GUIDE v2.5 17-Jun-2016 21:12:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AgarUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AgarUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



%%%%%%%%%%%%%%%%%%%%%%%
% Executes just before AgarUI is made visible.
% Loads the user file if it exists
% If there's no user file nothing is shown until user is added
% Pauses the UI until Agarstart is clicked

function AgarUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

if exist('Users.dat','file')
    handles.userinfo = table2cell(readtable('Users.dat'));
    handles.popup_selectuser.String = handles.userinfo(:,1);
    handles.nouserfile = 0;
    set(handles.box_ddataset,'visible','on')
    set(handles.box_viewchannels,'visible','on')
    set(handles.agarstart,'visible','on')
    handles.OSslash = handles.userinfo{1,2};
    handles.kymopath = handles.userinfo{1,3};
    handles.agarpath = strcat(handles.kymopath,'Agarolysis',handles.OSslash);
    handles.datapath = handles.kymopath;
else
    handles.nouserfile = 1;
    handles.datapath = ''; 
    handles.kymopath = '';
    handles.agarpath = '';
    handles.OSslash = '\';
end
handles.beampath = '';


% Choose default command line output for AgarUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AgarUI wait for user response (see UIRESUME)
uiwait(handles.figure1);



%%%%%%%%%%%%%%%%%%%%%%%
% Outputs from this function are returned to the command line.
% Stores the input values to the outputs of this function which is retured
% to AgarDefine

function varargout = AgarUI_OutputFcn(hObject, eventdata, handles) 

if handles.nouserfile
    msgbox('No user defined, add new user','Error','error')
end

% Get default command line output from handles structure

init.OSslash = handles.OSslash;
init.kymopath = handles.kymopath;
init.Agarpath = handles.agarpath;
init.datapath = handles.datapath;

init.PCimgname = handles.edit_PC.String;
init.CFPimgname = handles.edit_CFP.String;
init.YFPimgname = handles.edit_YFP.String;
init.RFPimgname = handles.edit_RFP.String;
init.meshfile = handles.edit_oufti.String;
init.meshpath = strcat(init.datapath,init.meshfile);

init.flresize = str2double(handles.edit_fl_resize.String);
init.fltrans = [str2double(handles.edit_fl_x.String),...
    str2double(handles.edit_fl_y.String)];
init.pcresize = str2double(handles.edit_pc_resize.String);
init.pctrans = [str2double(handles.edit_pc_x.String),...
    str2double(handles.edit_pc_y.String)];

init.CFPbeampath = strcat(handles.beampath,handles.edit_cfpbeam.String);
init.YFPbeampath = strcat(handles.beampath,handles.edit_yfpbeam.String);
init.RFPbeampath = strcat(handles.beampath,handles.edit_rfpbeam.String);

viewchanbool(1) = handles.checkbox_CFP.Value;
viewchanbool(2) = handles.checkbox_YFP.Value;
viewchanbool(3) = handles.checkbox_RFP.Value;

j = 1;
for i = 1:3; 
    if viewchanbool(i) == 1;
        init.viewchannels(j) = i;
        j = j+1;
    end
end

varargout{1} = init;
varargout{2} = handles.radio_ddata_yes.Value;
varargout{3} = handles.radio_dtrans_yes.Value;
varargout{4} = handles.radio_dbeam_yes.Value;

close(handles.figure1)



%%%%%%%%%%%%%%%%%%%%%%%
% Executes on button press in agarstart.
% Checks whether default files and beamshapes are selected (defined in
% AgarDefine). 
% If the default files are not selected, and the boxes are empty, errors
% will be prompted and the UI won't preceed. If there are no errors,
% uiresume will close the UI and continue the script. 

function agarstart_Callback(hObject, eventdata, handles)

ddata = logical(handles.radio_ddata_no.Value);
dbeam = logical(handles.radio_dbeam_no.Value);

if ddata
    file0 = isempty(handles.edit_folder.String)||strcmp(handles.edit_folder.String,'0');
    file1 = isempty(handles.edit_PC.String)||strcmp(handles.edit_PC.String,'0');
    file2 = isempty(handles.edit_CFP.String)||strcmp(handles.edit_CFP.String,'0');
    file3 = isempty(handles.edit_YFP.String)||strcmp(handles.edit_YFP.String,'0');
    file4 = isempty(handles.edit_RFP.String)||strcmp(handles.edit_RFP.String,'0');
    file5 = isempty(handles.edit_oufti.String)||strcmp(handles.edit_oufti.String,'0');

    file_error = file0 || file1 || file2 || file3 || file4 || file5;
    
    if dbeam
        file6 = isempty(handles.edit_cfpbeam.String)||strcmp(handles.edit_cfpbeam.String,'0');
        file7 = isempty(handles.edit_yfpbeam.String)||strcmp(handles.edit_yfpbeam.String,'0');
        file8 = isempty(handles.edit_rfpbeam.String)||strcmp(handles.edit_rfpbeam.String,'0');

        beam_error = file6 || file7 || file8;
    else
        beam_error = 0;
    end
else
    file_error = 0;
    beam_error = 0;
end

if file_error && beam_error
    msgbox('Files and Beamshapes missing','Error','error')
elseif file_error && ~beam_error
    msgbox('Files missing','Error','error')
elseif ~file_error && beam_error
    msgbox('Beamshapes missing','Error','error')
end

if ~(file_error || beam_error)
    uiresume(handles.figure1)
end


% --- Executes during object creation, after setting all properties.
function popup_selectuser_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject,handles);

% --- Executes on selection change in popup_selectuser.
function popup_selectuser_Callback(hObject, eventdata, handles)
selection = get(hObject,'Value');
handles.userinfo = table2cell(readtable('Users.dat'));
handles.OSslash = handles.userinfo{selection,2};
handles.kymopath = handles.userinfo{selection,3};
handles.agarpath = strcat(handles.kymopath,'Agarolysis',handles.OSslash);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New user UI


% --- Executes on button press in newuser.
function newuser_Callback(hObject, eventdata, handles)
set(handles.box_newuser,'visible','on')


function newuser_username_Callback(hObject, eventdata, handles)
% hObject    handle to newuser_username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of newuser_username as text
%        str2double(get(hObject,'String')) returns contents of newuser_username as a double


% --- Executes during object creation, after setting all properties.
function newuser_username_CreateFcn(hObject, eventdata, handles)
% hObject    handle to newuser_username (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in newuser_osbox.
function newuser_osbox_SelectionChangedFcn(hObject, eventdata, handles)
if get(handles.newuser_os_windows,'Value') == 1;
    handles.OSslash = '\';
elseif get(handles.newuser_os_mac,'Value') == 1;
    handles.OSslash = '/';
end
handles.agarpath = strcat(handles.kymopath,'Agarolysis',handles.OSslash);

guidata(hObject,handles)

% --- Executes on button press in newuser_kymopath.
function newuser_kymopath_Callback(hObject, eventdata, handles)

kymopath = uigetdir('','Select kymocode folder');
handles.kymopath = strcat(kymopath,handles.OSslash);
handles.agarpath = strcat(handles.kymopath,'Agarolysis',handles.OSslash);
handles.datapath = handles.kymopath;

guidata(hObject,handles)


% --- Executes on button press in newuser_done.
function newuser_done_Callback(hObject, eventdata, handles)
fname = ~isempty(get(handles.newuser_username,'String'));
fkymo = ~isempty(handles.kymopath);

if fname && fkymo
    kymopath = handles.kymopath;
    osslash = handles.OSslash;
    name = get(handles.newuser_username,'String');
    
    newinfo = {name,osslash,kymopath};
    
    if exist('Users.dat','file')
        userinfo = table2cell(readtable('Users.dat'));
        newinfo = [userinfo; newinfo];
    end
    
    writetable(cell2table(newinfo),'Users.dat')
    set(handles.box_newuser,'visible','off')
    set(handles.popup_selectuser,'String',newinfo(:,1));
    set(handles.box_ddataset,'visible','on')
    set(handles.box_viewchannels,'visible','on')
    set(handles.agarstart,'visible','on')
else
    msgbox('Missing information','Error','error') 
end


% --- Executes on button press in newuser_cancel.
function newuser_cancel_Callback(hObject, eventdata, handles)
set(handles.box_newuser,'visible','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default data or not

% --- Executes on button press in radio_ddata_yes.
function radio_ddata_yes_Callback(hObject, eventdata, handles)
set(handles.box_selectfiles,'visible','off');
set(handles.box_dtranslations,'visible','off');
set(handles.box_dbeamshapes,'visible','off');
set(handles.box_dtranslations,'visible','off');
set(handles.box_dbeamshapes,'visible','off');
set(handles.box_selectbeam,'visible','off');
set(handles.box_translations,'visible','off');
set(handles.save_preset,'visible','off');
set(handles.load_preset,'visible','off');
handles.radio_dtrans_yes.Value = 1;
handles.radio_dbeam_yes.Value = 1;


% --- Executes on button press in radio_ddata_no.
function radio_ddata_no_Callback(hObject, eventdata, handles)
set(handles.box_selectfiles,'visible','on');
set(handles.box_dtranslations,'visible','on');
set(handles.box_dbeamshapes,'visible','on');
set(handles.save_preset,'visible','on');
set(handles.load_preset,'visible','on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Presets

% --- Executes on button press in save_preset.
function save_preset_Callback(hObject, eventdata, handles)

data{1,1} = handles.edit_PC.String;
data{1,2} = handles.edit_CFP.String;
data{1,3} = handles.edit_YFP.String;
data{1,4} = handles.edit_RFP.String;
data{1,5} = handles.edit_oufti.String;
data{1,6} = handles.edit_fl_resize.String;
data{1,7} = handles.edit_fl_x.String;
data{1,8} = handles.edit_fl_y.String;
data{1,9} = handles.edit_pc_resize.String;
data{1,10} = handles.edit_pc_x.String;
data{1,11} = handles.edit_pc_y.String;
data{1,12} = handles.edit_cfpbeam.String;
data{1,13} = handles.edit_yfpbeam.String;
data{1,14} = handles.edit_rfpbeam.String;
data{1,15} = num2str(handles.checkbox_CFP.Value);
data{1,16} = num2str(handles.checkbox_YFP.Value);
data{1,17} = num2str(handles.checkbox_RFP.Value);
data{1,18} = handles.datapath;

filename = inputdlg('Set preset name','');
savename = strcat(filename{1},'.dat');
cd(handles.datapath)
writetable(cell2table(data),savename);
cd(handles.agarpath)


% --- Executes on button press in load_preset.
function load_preset_Callback(hObject, eventdata, handles)
[name,folder] = uigetfile('*.dat','Select Preset');
cd(folder)
data = table2cell(readtable(name));

set(handles.edit_PC,'String',data{1,1})
set(handles.edit_CFP,'String',data{1,2})
set(handles.edit_YFP,'String',data{1,3});
set(handles.edit_RFP,'String',data{1,4});
set(handles.edit_oufti,'String',data{1,5});
set(handles.edit_fl_resize,'String',data{1,6});
set(handles.edit_fl_x,'String',data{1,7});
set(handles.edit_fl_y,'String',data{1,8});
set(handles.edit_pc_resize,'String',data{1,9});
set(handles.edit_pc_x,'String',data{1,10});
set(handles.edit_pc_y,'String',data{1,11});
set(handles.edit_cfpbeam,'String',data{1,12});
set(handles.edit_yfpbeam,'String',data{1,13});
set(handles.edit_rfpbeam,'String',data{1,14});
set(handles.checkbox_CFP,'Value',data{1,15});
set(handles.checkbox_YFP,'Value',data{1,16});
set(handles.checkbox_RFP,'Value',data{1,17});
set(handles.edit_folder,'String',data{1,18});
handles.datapath = data{1,18};

cd(handles.agarpath)
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Files selection

% --- Executes on button press in select_folder.
function select_folder_Callback(hObject, eventdata, handles)
cd(handles.kymopath)
name = strcat(uigetdir('','Select data folder'),handles.OSslash);
cd(handles.agarpath)
set(handles.edit_folder,'String',name);
handles.datapath = name;
guidata(hObject,handles)

function select_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in select_PC.
function select_PC_Callback(hObject, eventdata, handles)
cd(handles.datapath)
[name,~] = uigetfile('*.tif','Select PC image');
cd(handles.agarpath)
set(handles.edit_PC,'String',name);
guidata(hObject,handles)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
cd(handles.datapath)
[name,~] = uigetfile('*.tif','Select CFP image');
cd(handles.agarpath)
set(handles.edit_CFP,'String',name);
guidata(hObject,handles)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
cd(handles.datapath)
[name,~] = uigetfile('*.tif','Select YFP image');
cd(handles.agarpath)
set(handles.edit_YFP,'String',name);
guidata(hObject,handles)

% --- Executes on button press in select_RFP.
function select_RFP_Callback(hObject, eventdata, handles)
cd(handles.datapath)
[name,~] = uigetfile('*.tif','Select RFP image');
cd(handles.agarpath)
set(handles.edit_RFP,'String',name);
guidata(hObject,handles)

% --- Executes on button press in select_oufti.
function select_oufti_Callback(hObject, eventdata, handles)
cd(handles.datapath)
[name,~] = uigetfile('*.mat','Select oufti output');
cd(handles.agarpath)
set(handles.edit_oufti,'String',name);
guidata(hObject,handles)




function edit_folder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_folder as text
%        str2double(get(hObject,'String')) returns contents of edit_folder as a double

% --- Executes during object creation, after setting all properties.
function edit_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.



function edit_PC_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_PC as text
%        str2double(get(hObject,'String')) returns contents of edit_PC as a double

% --- Executes during object creation, after setting all properties.
function edit_PC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CFP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_CFP as text
%        str2double(get(hObject,'String')) returns contents of edit_CFP as a double

% --- Executes during object creation, after setting all properties.
function edit_CFP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_YFP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_YFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_YFP as text
%        str2double(get(hObject,'String')) returns contents of edit_YFP as a double

% --- Executes during object creation, after setting all properties.
function edit_YFP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_YFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_RFP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_RFP as text
%        str2double(get(hObject,'String')) returns contents of edit_RFP as a double

% --- Executes during object creation, after setting all properties.
function edit_RFP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_oufti_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oufti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_oufti as text
%        str2double(get(hObject,'String')) returns contents of edit_oufti as a double

% --- Executes during object creation, after setting all properties.
function edit_oufti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oufti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default translations and default beamshapes

% --- Executes when selected object is changed in box_dtranslations.
function box_dtranslations_SelectionChangedFcn(hObject, eventdata, handles)
if get(handles.radio_dtrans_yes,'Value') == 1;
    set(handles.box_translations,'visible','off')
elseif get(handles.radio_dtrans_no,'Value') == 1;
    set(handles.box_translations,'visible','on')
end


% --- Executes when selected object is changed in box_dbeamshapes.
function box_dbeamshapes_SelectionChangedFcn(hObject, eventdata, handles)
if get(handles.radio_dbeam_yes,'Value') == 1;
    set(handles.box_selectbeam,'visible','off')
elseif get(handles.radio_dbeam_no,'Value') == 1;
    set(handles.box_selectbeam,'visible','on')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translations

function edit_fl_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fl_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_fl_x as text
%        str2double(get(hObject,'String')) returns contents of edit_fl_x as a double

% --- Executes during object creation, after setting all properties.
function edit_fl_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fl_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fl_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fl_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_fl_y as text
%        str2double(get(hObject,'String')) returns contents of edit_fl_y as a double

% --- Executes during object creation, after setting all properties.
function edit_fl_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fl_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fl_resize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fl_resize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_fl_resize as text
%        str2double(get(hObject,'String')) returns contents of edit_fl_resize as a double

% --- Executes during object creation, after setting all properties.
function edit_fl_resize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fl_resize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pc_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pc_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_pc_x as text
%        str2double(get(hObject,'String')) returns contents of edit_pc_x as a double

% --- Executes during object creation, after setting all properties.
function edit_pc_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pc_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pc_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pc_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_pc_y as text
%        str2double(get(hObject,'String')) returns contents of edit_pc_y as a double

% --- Executes during object creation, after setting all properties.
function edit_pc_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pc_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pc_resize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pc_resize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_pc_resize as text
%        str2double(get(hObject,'String')) returns contents of edit_pc_resize as a double

% --- Executes during object creation, after setting all properties.
function edit_pc_resize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pc_resize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% function edit_bf_x_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_bf_x (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Hints: get(hObject,'String') returns contents of edit_bf_x as text
% %        str2double(get(hObject,'String')) returns contents of edit_bf_x as a double
% 
% % --- Executes during object creation, after setting all properties.
% function edit_bf_x_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_bf_x (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% 
% function edit_bf_y_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_bf_y (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Hints: get(hObject,'String') returns contents of edit_bf_y as text
% %        str2double(get(hObject,'String')) returns contents of edit_bf_y as a double
% 
% % --- Executes during object creation, after setting all properties.
% function edit_bf_y_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_bf_y (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select beams

% --- Executes on button press in select_cfpbeam.
function select_cfpbeam_Callback(hObject, eventdata, handles)
if ~strcmp(handles.beampath,'');
    cd(handles.beampath)
end
[name,beampath] = uigetfile('*.tif','Select CFP beamshape');
cd(handles.agarpath)
set(handles.edit_cfpbeam,'String',name);
if ~strcmp(name,'0')
    handles.beampath = beampath;
end
guidata(hObject,handles)

% --- Executes on button press in select_yfpbeam
function select_yfpbeam_Callback(hObject, eventdata, handles)
if ~strcmp(handles.beampath,'');
    cd(handles.beampath)
end
[name,beampath] = uigetfile('*.tif','Select YFP beamshape');
cd(handles.agarpath)
set(handles.edit_yfpbeam,'String',name);
if ~strcmp(name,'0')
    handles.beampath = beampath;
end
guidata(hObject,handles)

% --- Executes on button press in select_rfpbeam.
function select_rfpbeam_Callback(hObject, eventdata, handles)
if ~strcmp(handles.beampath,'');
    cd(handles.beampath)
end
[name,beampath] = uigetfile('*.tif','Select RFP beamshape');
cd(handles.beampath)
set(handles.edit_rfpbeam,'String',name);
if ~strcmp(name,'0')
    handles.beampath = beampath;
end
guidata(hObject,handles)



function edit_rfpbeam_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rfpbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_rfpbeam as text
%        str2double(get(hObject,'String')) returns contents of edit_rfpbeam as a double

% --- Executes during object creation, after setting all properties.
function edit_rfpbeam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rfpbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yfpbeam_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yfpbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_yfpbeam as text
%        str2double(get(hObject,'String')) returns contents of edit_yfpbeam as a double

% --- Executes during object creation, after setting all properties.
function edit_yfpbeam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yfpbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cfpbeam_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cfpbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_cfpbeam as text
%        str2double(get(hObject,'String')) returns contents of edit_cfpbeam as a double

% --- Executes during object creation, after setting all properties.
function edit_cfpbeam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cfpbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View channels

% --- Executes on button press in checkbox_CFP.
function checkbox_CFP_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_CFP


% --- Executes on button press in checkbox_YFP.
function checkbox_YFP_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_YFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_YFP


% --- Executes on button press in checkbox_RFP.
function checkbox_RFP_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_RFP


% % --- Executes on button press in radio_data_yes.
% function radio_data_yes_Callback(hObject, eventdata, handles)
% 
% 
% % --- Executes on button press in radio_data_no.
% function radio_data_no_Callback(hObject, eventdata, handles)



% % --- Executes on button press in pushbutton16.
% function pushbutton16_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton16 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% 
% 
% function edit17_Callback(hObject, eventdata, handles)
% % hObject    handle to edit17 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Hints: get(hObject,'String') returns contents of edit17 as text
% %        str2double(get(hObject,'String')) returns contents of edit17 as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit17_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit17 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
