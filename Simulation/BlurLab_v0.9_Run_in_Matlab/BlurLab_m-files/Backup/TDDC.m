function varargout = TDDC(varargin)
% TDDC MATLAB code for TDDC.fig
%      TDDC, by itself, creates a new TDDC or raises the existing
%      singleton*.
%
%      H = TDDC returns the handle to a new TDDC or the handle to
%      the existing singleton*.
%
%      TDDC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TDDC.M with the given input arguments.
%
%      TDDC('Property','Value',...) creates a new TDDC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TDDC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TDDC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TDDC

% Last Modified by GUIDE v2.5 27-Sep-2016 15:54:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TDDC_OpeningFcn, ...
                   'gui_OutputFcn',  @TDDC_OutputFcn, ...
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


% --- Executes just before TDDC is made visible.
function TDDC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TDDC (see VARARGIN)


handles.plotx = [0, 1];
handles.ploty = [handles.slider_1.Value,handles.slider_2.Value];

Stringval1 = strcat(num2str(handles.slider_1.Value,'%0.f'),'%');
set(handles.text_s1,'String',Stringval1)
Stringval2 = strcat(num2str(handles.slider_2.Value,'%0.f'),'%');
set(handles.text_s2,'String',Stringval2)

axes(handles.axes1);
axis([0, 1, 0, 100]);
hold on
plot(handles.plotx,handles.ploty,'b--')
plot(handles.plotx,handles.ploty,'ro')
hold off

% Choose default command line output for TDDC
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

uiwait(handles.figure1)

% UIWAIT makes TDDC wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TDDC_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.plotx;
varargout{2} = handles.ploty;
varargout{3} = str2double(handles.edit_DC1.String);
varargout{4} = str2double(handles.edit_DC2.String);

if isnan(varargout{3}) || isnan(varargout{4})
    disp('Non numericals values entered for the diffusion constants')
end

close(handles.figure1)


% --- Executes on slider movement.
function slider_1_Callback(hObject, eventdata, handles)

handles.ploty(1) = handles.slider_1.Value;
Stringval = strcat(num2str(handles.slider_1.Value,'%0.f'),'%');
set(handles.text_s1,'String',Stringval)

axes(handles.axes1); 
cla(handles.axes1);

hold on
plot(handles.plotx,handles.ploty,'b--')
plot(handles.plotx,handles.ploty,'ro')
hold off

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_2_Callback(hObject, eventdata, handles)

handles.ploty(end) = handles.slider_2.Value;
Stringval = strcat(num2str(handles.slider_2.Value,'%0.f'),'%');
set(handles.text_s2,'String',Stringval)

axes(handles.axes1); 
cla(handles.axes1);

hold on
plot(handles.plotx,handles.ploty,'b--')
plot(handles.plotx,handles.ploty,'ro')
hold off

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_DC1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DC1 as text
%        str2double(get(hObject,'String')) returns contents of edit_DC1 as a double


% --- Executes during object creation, after setting all properties.
function edit_DC1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DC2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DC2 as text
%        str2double(get(hObject,'String')) returns contents of edit_DC2 as a double


% --- Executes during object creation, after setting all properties.
function edit_DC2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)

axes(handles.axes1); 
[x,y] = ginput(1);

if x > handles.plotx(end-1)
    
    points = numel(handles.plotx);
    handles.plotx = [handles.plotx(1:points-1),x,handles.plotx(end)];
    handles.ploty = [handles.ploty(1:points-1),y,handles.ploty(end)];

    cla(handles.axes1);

    hold on
    plot(handles.plotx,handles.ploty,'b--')
    plot(handles.plotx,handles.ploty,'ro')
    hold off
end

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in push_clear.
function push_clear_Callback(hObject, eventdata, handles)

handles.plotx = [handles.plotx(1), handles.plotx(end)];
handles.ploty = [handles.ploty(1), handles.ploty(end)];

axes(handles.axes1); 
cla(handles.axes1);

hold on
plot(handles.plotx,handles.ploty,'b--')
plot(handles.plotx,handles.ploty,'ro')
hold off

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in push_done.
function push_done_Callback(hObject, eventdata, handles)

uiresume(handles.figure1)
